#pragma once

#include <AMReX.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Print.H>

#include "Config.H"
#include "utils/HelperFunctions.H"
#include "utils/Constants.H"

#ifdef USE_IMAS
#include "io/ImasReader.H"
#include <AMReX_ParallelDescriptor.H>
#endif

/**
 * @brief Complete AMReX domain description for one simulation.
 *
 * `DomainSetup` stores both the integer index-space domain and the physical
 * coordinate-space `RealBox` used to construct an AMReX `Geometry`.
 */
struct DomainSetup
{
    /**
     * @brief Lower corner of the integer AMReX domain.
     */
    amrex::IntVect dom_lo;

    /**
     * @brief Upper corner of the integer AMReX domain.
     */
    amrex::IntVect dom_hi;

    /**
     * @brief Integer AMReX domain box.
     */
    amrex::Box domain;

    /**
     * @brief Physical coordinate bounds associated with the domain.
     */
    amrex::RealBox real_box;

    /**
     * @brief Periodicity flags for each coordinate direction.
     */
    std::array<int, AMREX_SPACEDIM> is_periodic;

    /**
     * @brief AMReX geometry built from the domain, real box, and periodicity.
     */
    amrex::Geometry geom;

    /**
     * @brief Construct a complete domain setup object.
     *
     * @param a_dom_lo Lower integer domain corner.
     * @param a_dom_hi Upper integer domain corner.
     * @param a_domain Integer domain box.
     * @param a_real_box Physical coordinate box.
     * @param a_is_periodic Periodicity flags.
     * @param a_geom AMReX geometry.
     */
    DomainSetup(const amrex::IntVect& a_dom_lo,
                const amrex::IntVect& a_dom_hi,
                const amrex::Box& a_domain,
                const amrex::RealBox& a_real_box,
                const std::array<int, AMREX_SPACEDIM>& a_is_periodic,
                const amrex::Geometry& a_geom)
        : dom_lo(a_dom_lo),
          dom_hi(a_dom_hi),
          domain(a_domain),
          real_box(a_real_box),
          is_periodic(a_is_periodic),
          geom(a_geom)
    {}
};

/**
 * @brief AMReX grid, distribution, and optional PML-shell setup.
 *
 * `GridSetup` stores the primary computational `BoxArray` and
 * `DistributionMapping`, as well as the exterior shell used by the PML.
 */
struct GridSetup
{
    /**
     * @brief BoxArray for the physical computational domain.
     */
    amrex::BoxArray ba;

    /**
     * @brief Distribution mapping for the physical computational domain.
     */
    amrex::DistributionMapping dm;

    /**
     * @brief BoxArray for the PML shell.
     */
    amrex::BoxArray ba_pml;

    /**
     * @brief Geometry for the PML shell.
     */
    amrex::Geometry geom_pml;

    /**
     * @brief Distribution mapping for the PML shell.
     */
    amrex::DistributionMapping dm_pml;
};

/**
 * @brief Build the AMReX domain from the runtime configuration.
 *
 * In manual-input mode, the domain is built from `input.n_cell`,
 * `input.prob_lo`, `input.prob_hi`, and `input.is_periodic`. In IMAS mode, the
 * domain is built from the precomputed `ImashSetup`.
 *
 * @param cfg Runtime configuration.
 * @param imas_setup Optional IMAS-derived setup, required when `cfg.use_imas`
 * is true.
 * @return Complete AMReX domain setup.
 */
inline DomainSetup make_domain_setup(const Config& cfg
#ifdef USE_IMAS
                                   , const ImashSetup* imas_setup = nullptr
#endif
)
{
    using namespace amrex;

    IntVect dom_lo(0, 0, 0);
    IntVect dom_hi;
    Box domain;
    RealBox real_box({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0});
    std::array<int, AMREX_SPACEDIM> is_per = {0, 0, 0};

#ifdef USE_IMAS
    if (cfg.use_imas)
    {
        /*
         * IMAS mode uses the cell counts and physical bounds already computed
         * by BuildImashSetup().
         */
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(imas_setup != nullptr,
            "make_domain_setup: cfg.use_imas=1 but imas_setup is null");

        dom_hi = IntVect((*imas_setup).n_cell[0] - 1,
                         (*imas_setup).n_cell[1] - 1,
                         (*imas_setup).n_cell[2] - 1);

        domain = Box(dom_lo, dom_hi);

        Real plo[AMREX_SPACEDIM] = {
            static_cast<Real>((*imas_setup).prob_lo[0]),
            static_cast<Real>((*imas_setup).prob_lo[1]),
            static_cast<Real>((*imas_setup).prob_lo[2])
        };

        Real phi[AMREX_SPACEDIM] = {
            static_cast<Real>((*imas_setup).prob_hi[0]),
            static_cast<Real>((*imas_setup).prob_hi[1]),
            static_cast<Real>((*imas_setup).prob_hi[2])
        };

        real_box = RealBox(plo, phi);

        is_per = {
            (*imas_setup).is_periodic[0],
            (*imas_setup).is_periodic[1],
            (*imas_setup).is_periodic[2]
        };

        Geometry geom(domain, &real_box, CoordSys::cartesian, is_per.data());

        return DomainSetup(dom_lo, dom_hi, domain, real_box, is_per, geom);
    }
#endif

    /*
     * Manual mode starts from the configured cell counts. Directions with more
     * than one cell are multiplied by input_scale; one-cell directions remain
     * collapsed.
     */
    int szx = cfg.input_n_cell[0];
    int szy = cfg.input_n_cell[1];
    int szz = cfg.input_n_cell[2];

    if (szx > 1) { szx *= cfg.input_scale; }
    if (szy > 1) { szy *= cfg.input_scale; }
    if (szz > 1) { szz *= cfg.input_scale; }

    dom_hi = IntVect(--szx, --szy, --szz);

    domain = Box(dom_lo, dom_hi);

    Real plo[AMREX_SPACEDIM] = {
        static_cast<Real>(cfg.input_prob_lo[0]),
        static_cast<Real>(cfg.input_prob_lo[1]),
        static_cast<Real>(cfg.input_prob_lo[2])
    };

    Real phi[AMREX_SPACEDIM] = {
        static_cast<Real>(cfg.input_prob_hi[0]),
        static_cast<Real>(cfg.input_prob_hi[1]),
        static_cast<Real>(cfg.input_prob_hi[2])
    };

    real_box = RealBox(plo, phi);

    is_per = {
        cfg.input_is_periodic[0],
        cfg.input_is_periodic[1],
        cfg.input_is_periodic[2]
    };

    Geometry geom(domain, &real_box, CoordSys::cartesian, is_per.data());

    return DomainSetup(dom_lo, dom_hi, domain, real_box, is_per, geom);
}

/**
 * @brief Build the physical and PML AMReX grid layouts.
 *
 * The physical domain is boxed and optionally split using `max_grid_size`. The
 * PML shell is then generated around the physical domain using
 * `outer_shell_generation()`.
 *
 * @param cfg Runtime configuration.
 * @param dom Domain setup.
 * @return Grid setup for the physical domain and PML shell.
 */
inline GridSetup make_grid_setup(const Config& cfg,
                                 const DomainSetup& dom)
{
    using namespace amrex;

    BoxArray ba(dom.domain);

    if (cfg.max_grid_size > 0) {
        ba.maxSize(cfg.max_grid_size);
    }

    DistributionMapping dm(ba);

    BoxArray ba_pml;
    Geometry geom_pml;
    DistributionMapping dm_pml;

    outer_shell_generation(ba,
                           dom.geom,
                           ba_pml,
                           geom_pml,
                           dm_pml,
                           cfg.pml_thickness);

    return GridSetup{ba, dm, ba_pml, geom_pml, dm_pml};
}

/**
 * @brief Compute the maximum stable electromagnetic timestep.
 *
 * The estimate uses the Cartesian FDTD CFL condition. When `dimensional` is
 * true, the physical speed of light is included. When false, normalized units
 * with `c = 1` are assumed.
 *
 * @param geom AMReX geometry defining the cell sizes.
 * @param dimensional If true, include `c_0` in the CFL estimate.
 * @return Maximum stable timestep.
 */
inline amrex::Real compute_dt_max(const amrex::Geometry& geom,
                                  bool dimensional)
{
    using namespace amrex;

    const auto dx = geom.CellSizeArray();

    Real denom = std::sqrt(
        1.0 / (dx[0] * dx[0]) +
        1.0 / (dx[1] * dx[1]) +
        1.0 / (dx[2] * dx[2])
    );

    return dimensional ? 1.0 / (c_0 * denom)
                       : 1.0 / denom;
}

/**
 * @brief Print characteristic scales, runtime options, domain data, and
 * timestep information.
 *
 * If unit normalization is active, dimensional values are reconstructed for
 * human-readable output.
 *
 * @param cfg Runtime configuration.
 * @param dom Domain setup.
 * @param dt_max Maximum stable timestep.
 * @param dt Selected timestep.
 */
inline void print_config_summary(const Config& cfg,
                                 const DomainSetup& dom,
                                 amrex::Real dt_max,
                                 amrex::Real dt)
{
    using namespace amrex;

    Real lg = 1.0;
    Real tm = 1.0;

    if (cfg.normalize_units) {
        lg = L_0;
        tm = t_0;
    }

    Print() << "\n====== Characteristic scales ======\n";
    Print() << "distance (m)           = " << L_0 << "\n";
    Print() << "B field (T)            = " << B_0 << "\n";
    Print() << "Density (1/m^3)        = " << n_0 << "\n\n";

    Print() << "\n========== Derived scales ==========\n";
    Print() << "time (s)               = " << t_0 << "\n";
    Print() << "Temperature (V)        = " << T_0 << "\n";
    Print() << "Charge density (C/m^3) = " << rho_0 << "\n";
    Print() << "H field (A/m)          = " << H_0 << "\n";
    Print() << "E field (V/m)          = " << E_0 << "\n";
    Print() << "Current J (A/m^2)      = " << J_0 << "\n";
    Print() << "Collision F (Hz)       = " << nu_0 << "\n";
    Print() << "Charge (C)             = " << q_0 << "\n";
    Print() << "Mass (kg)              = " << m_0 << "\n\n";

    Print() << "\n======= Configuration summary =======\n";
    Print() << "output_dir             = " << cfg.output_dir << "\n";
    Print() << "max_grid_size          = " << cfg.max_grid_size << "\n";
    Print() << "nghost                 = " << cfg.nghost << "\n";
    Print() << "pml_thickness          = " << cfg.pml_thickness << "\n";
    Print() << "cfl                    = " << cfg.cfl << "\n";
    Print() << "final_time             = " << cfg.final_time * tm << "\n";
    Print() << "max_step               = " << cfg.max_step << "\n";

#ifdef USE_IMAS
    if (cfg.use_imas)
    {
        /*
         * IMAS mode reports the domain actually generated from the IMAS setup.
         */
        Print() << "domain cells           = "
                << dom.domain.length(0) << " x "
                << dom.domain.length(1) << " x "
                << dom.domain.length(2) << "\n";

        Print() << "prob_lo                = "
                << dom.real_box.lo(0) * lg << " "
                << dom.real_box.lo(1) * lg << " "
                << dom.real_box.lo(2) * lg << "\n";

        Print() << "prob_hi                = "
                << dom.real_box.hi(0) * lg << " "
                << dom.real_box.hi(1) * lg << " "
                << dom.real_box.hi(2) * lg << "\n";

        Print() << "is_periodic            = "
                << dom.is_periodic[0] << " "
                << dom.is_periodic[1] << " "
                << dom.is_periodic[2] << "\n";

        Print() << "n_cell                 = "
                << cfg.input_n_cell[0] << " "
                << cfg.input_n_cell[1] << " "
                << cfg.input_n_cell[2] << "\n";
    }
    else
#endif
    {
        /*
         * Manual mode reports the user-provided physical box and periodicity.
         */
        Print() << "domain cells           = "
                << dom.domain.length(0) << " x "
                << dom.domain.length(1) << " x "
                << dom.domain.length(2) << "\n";

        Print() << "prob_lo                = "
                << cfg.input_prob_lo[0] * lg << " "
                << cfg.input_prob_lo[1] * lg << " "
                << cfg.input_prob_lo[2] * lg << "\n";

        Print() << "prob_hi                = "
                << cfg.input_prob_hi[0] * lg << " "
                << cfg.input_prob_hi[1] * lg << " "
                << cfg.input_prob_hi[2] * lg << "\n";

        Print() << "is_periodic            = "
                << cfg.input_is_periodic[0] << " "
                << cfg.input_is_periodic[1] << " "
                << cfg.input_is_periodic[2] << "\n";
    }

    Print() << "dt_max / t_0           = " << dt_max << "\n";
    Print() << "dt / t_0               = " << dt << " (< 10 for PML)\n";
    Print() << "================================\n\n";
}

#ifdef USE_IMAS

/**
 * @brief Convert an IMAS setup from physical units to normalized solver units.
 *
 * @param setup IMAS setup to scale in place.
 */
inline void remove_dimension(ImashSetup& setup)
{
    setup.prob_lo[0] /= L_0;
    setup.prob_lo[1] /= L_0;
    setup.prob_lo[2] /= L_0;

    setup.prob_hi[0] /= L_0;
    setup.prob_hi[1] /= L_0;
    setup.prob_hi[2] /= L_0;

    setup.dx /= L_0;
    setup.dy /= L_0;
    setup.dz /= L_0;

    setup.dt_max /= t_0;
    setup.dt     /= t_0;
}

/**
 * @brief Restore an IMAS setup from normalized solver units to physical units.
 *
 * @param setup IMAS setup to scale in place.
 */
inline void add_dimension(ImashSetup& setup)
{
    setup.prob_lo[0] *= L_0;
    setup.prob_lo[1] *= L_0;
    setup.prob_lo[2] *= L_0;

    setup.prob_hi[0] *= L_0;
    setup.prob_hi[1] *= L_0;
    setup.prob_hi[2] *= L_0;

    setup.dx *= L_0;
    setup.dy *= L_0;
    setup.dz *= L_0;

    setup.dt_max *= t_0;
    setup.dt     *= t_0;
}

/**
 * @brief Broadcast one integer from the IO processor to all MPI ranks.
 *
 * @param x Integer value to broadcast in place.
 */
inline void bcast_int(int& x)
{
    amrex::ParallelDescriptor::Bcast(&x, 1, 0);
}

/**
 * @brief Broadcast one double from the IO processor to all MPI ranks.
 *
 * @param x Double value to broadcast in place.
 */
inline void bcast_double(double& x)
{
    amrex::ParallelDescriptor::Bcast(&x, 1, 0);
}

/**
 * @brief Broadcast a vector of doubles from the IO processor to all MPI ranks.
 *
 * Non-IO ranks resize their vector before receiving the data.
 *
 * @param v Vector to broadcast in place.
 */
inline void bcast_vector_double(std::vector<double>& v)
{
    int n = static_cast<int>(v.size());

    bcast_int(n);

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        v.resize(n);
    }

    if (n > 0) {
        amrex::ParallelDescriptor::Bcast(v.data(), n, 0);
    }
}

/**
 * @brief Broadcast an axisymmetric R-Z grid.
 *
 * @param g Grid to broadcast in place.
 */
inline void BcastAxisymmetricGrid2D(AxisymmetricGrid2D& g)
{
    bcast_int(g.nR);
    bcast_int(g.nZ);

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        g.Rvals.resize(g.nR);
        g.Zvals.resize(g.nZ);
    }

    bcast_vector_double(g.Rvals);
    bcast_vector_double(g.Zvals);
}

/**
 * @brief Broadcast a scalar field tabulated on an R-Z grid.
 *
 * @param f Scalar field to broadcast in place.
 */
inline void BcastScalarField2D(ScalarField2D& f)
{
    bcast_int(f.nR);
    bcast_int(f.nZ);

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        f.values.resize(f.nR * f.nZ);
    }

    bcast_vector_double(f.values);
}

/**
 * @brief Broadcast cylindrical magnetic-field data on an R-Z grid.
 *
 * @param B Vector magnetic-field data to broadcast in place.
 */
inline void BcastVectorFieldTokamak2D(VectorFieldTokamak2D& B)
{
    bcast_int(B.nR);
    bcast_int(B.nZ);

    if (!amrex::ParallelDescriptor::IOProcessor()) {
        B.BR.resize(B.nR * B.nZ);
        B.BZ.resize(B.nR * B.nZ);
        B.Bphi.resize(B.nR * B.nZ);
    }

    bcast_vector_double(B.BR);
    bcast_vector_double(B.BZ);
    bcast_vector_double(B.Bphi);
}

/**
 * @brief Broadcast all plasma profile data needed by FDTD initialization.
 *
 * @param p Plasma data to broadcast in place.
 */
inline void BcastPlasmaData(PlasmaData& p)
{
    BcastAxisymmetricGrid2D(p.grid);

    BcastScalarField2D(p.rho);

    BcastScalarField2D(p.ne_core);
    BcastScalarField2D(p.Te_core);

    BcastScalarField2D(p.ne_edge);
    BcastScalarField2D(p.Te_edge);

    BcastScalarField2D(p.ne);
    BcastScalarField2D(p.Te);

    BcastVectorFieldTokamak2D(p.B);

    BcastScalarField2D(p.core_mask);
    BcastScalarField2D(p.edge_mask);
    BcastScalarField2D(p.source_mask);
}

/**
 * @brief Broadcast an IMAS setup from the IO processor to all MPI ranks.
 *
 * Wall geometry is intentionally not broadcast because the FDTD plasma fill
 * only requires plasma profiles, magnetic-axis metadata, grid bounds, and
 * timestep/domain metadata.
 *
 * @param s IMAS setup to broadcast in place.
 */
inline void BcastImashSetup(ImashSetup& s)
{
    BcastPlasmaData(s.plasma);

    bcast_double(s.axis.R);
    bcast_double(s.axis.Z);

    bcast_int(s.plasma_nR);
    bcast_int(s.plasma_nZ);

    bcast_double(s.magnetic_axis_R);
    bcast_double(s.magnetic_axis_Z);

    bcast_int(s.n_wall_contours);
    bcast_int(s.n_wall_meshes);

    bcast_int(s.is_3d);

    for (int d = 0; d < 3; ++d) {
        bcast_double(s.prob_lo[d]);
        bcast_double(s.prob_hi[d]);
        bcast_int(s.n_cell[d]);
        bcast_int(s.is_periodic[d]);
    }

    bcast_double(s.dx);
    bcast_double(s.dy);
    bcast_double(s.dz);

    bcast_double(s.dt_max);
    bcast_double(s.dt);
}

#endif