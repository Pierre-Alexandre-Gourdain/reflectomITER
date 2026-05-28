#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>

#include <iomanip>
#include <memory>
#include <vector>
#include <cstdlib>
#include <string>

#include "io/Config.H"
#include "io/ConfigHelpers.H"

#ifdef USE_IMAS
#include "io/ImasReader.H"
#endif

#include "fdtd/EHJSolver.H"
#include "fdtd/PMLSolver.H"
#include "sources/ADEPlasma.H"
#include "fdtd/DivCleaning.H"
#include "sources/PointSource.H"
#include "utils/HelperFunctions.H"

using namespace amrex;
using namespace amrex::ParallelDescriptor;

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    double start_wall_time = ParallelDescriptor::second();
    double time_scale = 1.0;

    {
        /*
         * Read and validate the runtime configuration before any grid or solver
         * objects are constructed. Validation catches inconsistent input early,
         * before AMReX data structures are allocated.
         */
        Config cfg;
        cfg.read();
        cfg.validate();

#ifdef USE_IMAS
        /*
         * In IMAS mode, rank zero reads the IMAS-backed plasma/equilibrium data
         * and constructs the derived domain metadata. The result is then
         * broadcast so every MPI rank can build the same AMReX domain and sample
         * plasma data consistently.
         */
        std::unique_ptr<ImashSetup> imas_setup;

        if (cfg.use_imas) {
            imas_setup = std::make_unique<ImashSetup>();

            if (amrex::ParallelDescriptor::IOProcessor()) {
                *imas_setup = BuildImashSetup(cfg);
            }

            BcastImashSetup(*imas_setup);
        }
#endif

#ifdef USE_IMAS
        /*
         * Unit normalization is applied after the IMAS setup is built so that
         * IMAS extraction uses physical units, while the solver can still run in
         * normalized units.
         */
        if (cfg.use_imas) {
            if (cfg.normalize_units) {
                cfg.remove_dimension();
                remove_dimension(*imas_setup);
            }
        } else
#endif
        {
            /*
             * Manual-input runs can be normalized directly after configuration
             * validation.
             */
            if (cfg.normalize_units) {
                cfg.remove_dimension();
            }
        }

        /*
         * Build the AMReX domain and the associated physical/PML BoxArrays. The
         * domain may come either from manual input or from the IMAS-derived
         * setup.
         */
        DomainSetup dom =
#ifdef USE_IMAS
            make_domain_setup(cfg, imas_setup.get());
#else
            make_domain_setup(cfg);
#endif

        GridSetup grid = make_grid_setup(cfg, dom);

        /*
         * The timestep is computed from the physical speed of light unless the
         * configuration has already been normalized to code units.
         */
        Real dt_max = compute_dt_max(dom.geom, !cfg.normalize_units);
        Real dt = cfg.cfl * dt_max;

        Real time = 0.0;
        Real final_time = cfg.final_time;

        print_config_summary(cfg, dom, dt_max, dt);

        /*
         * The primary solver owns the physical-domain electromagnetic fields
         * and the current density used by the E update.
         */
        EHJSolver solver(grid.ba, dom.geom, grid.dm, cfg.nghost);
        solver.output_dir = cfg.output_dir;

        /*
         * Optional physics and numerical-control components are allocated only
         * when enabled by the input file. unique_ptr keeps ownership local and
         * makes the main loop checks explicit.
         */
        std::unique_ptr<ADEPlasma>   plasma;
        std::unique_ptr<ADEPlasma>   pml_plasma;
        std::unique_ptr<PMLSolver>   pml_solver;
        std::unique_ptr<DivCleaning> E_cleaner;
        std::unique_ptr<DivCleaning> H_cleaner;

        if (cfg.use_pml) {
            pml_solver = std::make_unique<PMLSolver>(
                solver,
                grid.ba_pml,
                grid.geom_pml,
                grid.dm_pml,
                cfg.pml_thickness,
                std::abs(cfg.pml_sigma)
            );

            pml_solver->output_dir = cfg.output_dir;
        }

        if (cfg.use_e_cleaner) {
            E_cleaner = std::make_unique<DivCleaning>(solver);
            E_cleaner->output_dir = cfg.output_dir;
        }

        if (cfg.use_plasma) {
            /*
             * The plasma model accumulates current into the attached solver.
             * When a PML exists, a second plasma object is attached to the PML
             * solver so source/current handling remains consistent in the
             * absorbing layer.
             */
            plasma = std::make_unique<ADEPlasma>(solver);
            plasma->output_dir = cfg.output_dir;

            pml_plasma = std::make_unique<ADEPlasma>(*pml_solver);
            pml_plasma->output_dir = cfg.output_dir;

#ifdef USE_IMAS
            /*
             * IMAS plasma data are sampled onto the physical and PML grids.
             */
            if (cfg.use_imas && plasma) {
                FillPlasmaFromImas(*imas_setup, *plasma, dom.geom);
            }

            if (cfg.use_imas && pml_plasma) {
                FillPlasmaFromImas(*imas_setup, *pml_plasma, dom.geom);
            }
#endif
        }

        if (cfg.use_h_cleaner) {
            H_cleaner = std::make_unique<DivCleaning>(solver, true);
            H_cleaner->output_dir = cfg.output_dir;
        }

        if (cfg.normalize_units)
        {
            /*
             * The solver runs in normalized units, but plot output can still be
             * written in dimensional units by enabling each output manager's
             * unit-scaling path.
             */
            time_scale *= t_0;

            solver.remove_dimension();

            if (plasma)     plasma->remove_dimension();
            if (pml_plasma) pml_plasma->remove_dimension();
            if (pml_solver) pml_solver->remove_dimension();
            if (E_cleaner)  E_cleaner->remove_dimension();
            if (H_cleaner)  H_cleaner->remove_dimension();
        }

        /*
         * Convert enabled source configuration blocks into source objects. The
         * source constructors receive normalized or dimensional quantities
         * consistently with the rest of the configuration.
         */
        std::vector<std::unique_ptr<PointSource>> sources;
        sources.reserve(cfg.sources.size());

        auto parse_source_type = [](const std::string& s) -> Source::SourceType
        {
            if (s == "J" || s == "j") return Source::SourceType::J;
            if (s == "E" || s == "e") return Source::SourceType::E;
            if (s == "H" || s == "h") return Source::SourceType::H;

            amrex::Abort("Invalid source type: '" + s + "'. Allowed: J, E, H");
            return Source::SourceType::J;
        };

        auto parse_profile_type = [](const std::string& s) -> Source::ProfileType
        {
            if (s == "gaussian" || s == "Gaussian") {
                return Source::ProfileType::Gaussian;
            }

            if (s == "supergaussian" ||
                s == "SuperGaussian" ||
                s == "super_gaussian") {
                return Source::ProfileType::SuperGaussian;
            }

            if (s == "heaviside" || s == "Heaviside") {
                return Source::ProfileType::Heaviside;
            }

            amrex::Abort("Invalid source profile: '" + s +
                         "'. Allowed: gaussian, supergaussian, heaviside");
            return Source::ProfileType::Gaussian;
        };

        for (const auto& scfg : cfg.sources)
        {
            if (!scfg.enabled) continue;

            amrex::Vector<amrex::Real> direction(AMREX_SPACEDIM);
            amrex::Vector<amrex::Real> shape(AMREX_SPACEDIM);

            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                direction[d] = scfg.direction[d];
                shape[d] = scfg.shape[d];
            }

            sources.push_back(std::make_unique<PointSource>(
                solver,
                parse_source_type(scfg.type),
                scfg.position[0],
                scfg.position[1],
                scfg.position[2],
                scfg.amplitude,
                scfg.frequency,
                scfg.phase,
                scfg.width,
                parse_profile_type(scfg.profile),
                scfg.supergaussian_order,
                direction,
                dom.geom,
                shape,
                scfg.t_on,
                scfg.t_off,
                scfg.t_rise
            ));
        }

        /*
         * If both PML and plasma are active, use the PML scalar profile to set
         * the conductivity mask for the plasma object attached to the PML.
         */
        if (pml_solver && pml_plasma) {
            pml_plasma->fill_sigma(*pml_solver, cfg.plasma_sigma);
        }

        /*
         * Main explicit time loop. The order is:
         *
         *   1. reset current,
         *   2. accumulate plasma current,
         *   3. apply configured sources,
         *   4. exchange PML current,
         *   5. update E and optional E cleaning,
         *   6. update H and optional H cleaning,
         *   7. write output when requested.
         */
        int step = 0;
        int nout = 0;

        double print_time = final_time / cfg.number_of_outputs;
        double lap_time = 0.0;
        double step_time = ParallelDescriptor::second();

        while (time <= final_time &&
               (cfg.max_step < 0 || step <= cfg.max_step))
        {
            /*
             * Current density is accumulated from plasma and source terms each
             * step, so it is reset before the current-producing components run.
             */
            solver.j().setVal(0.0);

            if (pml_solver) {
                pml_solver->j().setVal(0.0);
            }

            if (plasma) {
                plasma->update(dt);
            }

            if (pml_plasma) {
                pml_plasma->update(dt);
            }

            for (auto& src : sources) {
                src->update(time);
            }

            if (pml_solver) {
                pml_solver->FullBoundaryExchangeJ();
            }

            /*
             * Electric update. The PML auxiliary field is advanced before the
             * PML electric update so the correction corresponds to the current
             * timestep.
             */
            if (pml_solver) {
                pml_solver->updatePhiE(dt);
            }

            solver.updateE(dt);

            if (pml_solver) {
                pml_solver->updateE(dt);
            }

            if (E_cleaner) {
                E_cleaner->updatePhi(dt);
                E_cleaner->correctE(dt);
            }

            /*
             * Magnetic update and optional magnetic divergence cleaning.
             */
            if (pml_solver) {
                pml_solver->updatePhiH(dt);
            }

            solver.updateH(dt);

            if (pml_solver) {
                pml_solver->updateH(dt);
            }

            if (H_cleaner) {
                H_cleaner->updatePhi(dt);
                H_cleaner->correctH(dt);
            }

            /*
             * Write output at the requested cadence and always at the initial
             * state. The output time is converted back to dimensional units when
             * unit normalization is active.
             */
            if (lap_time > print_time || time == 0) {
                Print() << "Output number " << nout
                        << ", step " << step
                        << " completed, time = " << time * time_scale
                        << ", dt = " << dt * time_scale
                        << ", Wall time = "
                        << ParallelDescriptor::second() - step_time
                        << " s\n";

                solver.writeAMReXPlotfile(nout, time * time_scale);

                if (pml_solver) {
                    pml_solver->writeAMReXPlotfile(nout, time * time_scale);
                }

                lap_time = 0.0;
                nout++;
                step_time = ParallelDescriptor::second();
            }

            time += dt;
            lap_time += dt;
            ++step;
        }

        /*
         * Copy the input file into the output directory for reproducibility.
         * This assumes argv[1] is the input file path used to launch the run.
         */
		 
		if (argc > 1) {
			std::string str = std::string("cp \"")
							+ argv[1]
							+ "\" \""
							+ cfg.output_dir
							+ "/outputs/"
							+ "\"";
			int ret = ::system(str.c_str());
			(void) ret;
		}
    }

    double stop_wall_time = ParallelDescriptor::second();

    Print() << "Total sim time :"
            << stop_wall_time - start_wall_time
            << " s\n";

    amrex::Finalize();

    return 0;
}