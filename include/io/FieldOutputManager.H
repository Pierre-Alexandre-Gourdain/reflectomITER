#ifndef FIELDOUTPUT_MANAGER_H_
#define FIELDOUTPUT_MANAGER_H_

#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>
#include <AMReX_BoxArray.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Vector.H>
#include <AMReX_RealBox.H>

#include <filesystem>
#include <fstream>
#include <string>
#include <iomanip>
#include <set>
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>

#include "utils/Constants.H"
#include "utils/HelperFunctions.H"

namespace fs = std::filesystem;

/**
 * @brief Manages field registration and AMReX plotfile output.
 *
 * `FieldOutputManager` provides a common output interface for solver and
 * diagnostic classes. Classes can register one or more AMReX `MultiFab` fields
 * with names, geometries, component ranges, unit scales, and level indices.
 * Registered fields can then be written as AMReX multilevel plotfiles.
 *
 * The class also supports optional unit scaling. When unit scaling is enabled,
 * field values and grid coordinates are multiplied by the scales supplied at
 * registration time before output.
 *
 * The destructor writes `.visit` file lists on MPI rank zero so time series of
 * AMReX plotfiles can be opened conveniently in visualization tools such as
 * VisIt.
 */
class FieldOutputManager
{
public:
    /**
     * @brief Directory where output data are written.
     *
     * The output routines create subdirectories under this path.
     */
    std::string output_dir;

    /**
     * @brief Enable or disable unit scaling during output.
     *
     * @param v If true, registered field and grid scales are applied when
     * writing output.
     */
    void set_apply_unit_scaling(bool v) { apply_unit_scaling = v; }

    /**
     * @brief Return whether unit scaling is currently enabled.
     *
     * @return True if output scaling is enabled.
     */
    bool uses_unit_scaling() const { return apply_unit_scaling; }

    /**
     * @brief Enable output scaling to dimensionless units.
     *
     * Derived classes may override this method to scale internal diagnostic
     * fields in addition to enabling the manager-level output scaling flag.
     */
    virtual void remove_dimension() { apply_unit_scaling = true; }

    /**
     * @brief Disable output scaling and restore dimensional output.
     *
     * Derived classes may override this method to restore internal diagnostic
     * fields in addition to disabling the manager-level output scaling flag.
     */
    virtual void add_dimension() { apply_unit_scaling = false; }

public:
    /**
     * @brief Metadata describing one registered output field.
     *
     * A `FieldData` entry stores the field name, the source `MultiFab`, the
     * geometry used to write it, the associated `BoxArray`, the component range
     * to output, normalization factors, and the AMReX refinement level.
     */
    struct FieldData
    {
        /**
         * @brief Base name used for the output variable.
         */
        std::string name;

        /**
         * @brief Pointer to the field data.
         *
         * The manager does not own this `MultiFab`; it must remain valid while
         * output is requested.
         */
        amrex::MultiFab* mf;

        /**
         * @brief Geometry associated with the field.
         *
         * The manager does not own this geometry; it must remain valid while
         * output is requested.
         */
        const amrex::Geometry* geom;

        /**
         * @brief BoxArray associated with the field.
         */
        amrex::BoxArray ba;

        /**
         * @brief First component included in output.
         */
        int i_min;

        /**
         * @brief One-past-last component included in output.
         */
        int i_max;

        /**
         * @brief AMReX refinement level associated with this field.
         */
        int level;

        /**
         * @brief Multiplicative scale applied to field values when unit scaling
         * is enabled.
         */
        amrex::Real scale;

        /**
         * @brief Multiplicative scale applied to grid coordinates when unit
         * scaling is enabled.
         */
        amrex::Real grid_scale;

        /**
         * @brief Register all components of a field for output.
         *
         * @param n Output variable base name.
         * @param m Source field.
         * @param g Geometry associated with the source field.
         * @param b BoxArray associated with the source field.
         * @param sc Field-value scale used when output scaling is enabled.
         * @param gsc Grid-coordinate scale used when output scaling is enabled.
         * @param l AMReX refinement level.
         */
        FieldData(const std::string& n,
                  amrex::MultiFab* m,
                  const amrex::Geometry* g,
                  const amrex::BoxArray& b,
                  amrex::Real sc = 1.,
                  amrex::Real gsc = 1.,
                  int l = 0)
            : name(n),
              mf(m),
              geom(g),
              ba(b),
              i_min(0),
              i_max(m->nComp()),
              level(l),
              scale(sc),
              grid_scale(gsc)
        {}

        /**
         * @brief Register a component subset of a field for output.
         *
         * @param n Output variable base name.
         * @param m Source field.
         * @param g Geometry associated with the source field.
         * @param b BoxArray associated with the source field.
         * @param _i_min First component included in output.
         * @param _i_max One-past-last component included in output.
         * @param sc Field-value scale used when output scaling is enabled.
         * @param gsc Grid-coordinate scale used when output scaling is enabled.
         * @param l AMReX refinement level.
         */
        FieldData(const std::string& n,
                  amrex::MultiFab* m,
                  const amrex::Geometry* g,
                  const amrex::BoxArray& b,
                  int _i_min,
                  int _i_max,
                  amrex::Real sc = 1.,
                  amrex::Real gsc = 1.,
                  int l = 0)
            : name(n),
              mf(m),
              geom(g),
              ba(b),
              i_min(_i_min),
              i_max(_i_max),
              level(l),
              scale(sc),
              grid_scale(gsc)
        {}
    };

    /**
     * @brief Destroy the manager and write `.visit` time-series files.
     *
     * Rank zero scans the output directory and writes one `.visit` file per
     * output series. This makes AMReX plotfile sequences easier to open in
     * visualization software.
     */
    virtual ~FieldOutputManager()
    {
        int rank = amrex::ParallelDescriptor::MyProc();

        if (!rank) {
            WriteAMReXMasterHeaderList();
        }
    }

    /**
     * @brief Set the output database/series name for this manager.
     *
     * The name is appended under `outputs/` when plotfiles are written.
     *
     * @param database_name_ Output series name.
     */
    void set_name(const std::string& database_name_)
    {
        database_name = '/' + database_name_;
    }

    /**
     * @brief Register all components of a `MultiFab` for output.
     *
     * @param name Output variable base name.
     * @param mf Field to output.
     * @param geom Geometry associated with the field.
     * @param ba BoxArray associated with the field.
     * @param scale Field-value scale used when output scaling is enabled.
     * @param grid_scale Grid-coordinate scale used when output scaling is
     * enabled.
     * @param level AMReX refinement level.
     */
    void addField(const std::string& name,
                  amrex::MultiFab& mf,
                  const amrex::Geometry& geom,
                  const amrex::BoxArray& ba,
                  amrex::Real scale = 1.,
                  amrex::Real grid_scale = 1.,
                  int level = 0)
    {
        fields.emplace_back(name, &mf, &geom, ba, scale, grid_scale, level);
    }

    /**
     * @brief Register a component subset of a `MultiFab` for output.
     *
     * @param name Output variable base name.
     * @param mf Field to output.
     * @param geom Geometry associated with the field.
     * @param ba BoxArray associated with the field.
     * @param i_min First component included in output.
     * @param i_max One-past-last component included in output.
     * @param scale Field-value scale used when output scaling is enabled.
     * @param grid_scale Grid-coordinate scale used when output scaling is
     * enabled.
     * @param level AMReX refinement level.
     */
    void addField(const std::string& name,
                  amrex::MultiFab& mf,
                  const amrex::Geometry& geom,
                  const amrex::BoxArray& ba,
                  int i_min,
                  int i_max,
                  amrex::Real scale,
                  amrex::Real grid_scale,
                  int level)
    {
        fields.emplace_back(name, &mf, &geom, ba, i_min, i_max,
                            scale, grid_scale, level);
    }

    /**
     * @brief Write all registered fields as an AMReX multilevel plotfile.
     *
     * Registered fields are flattened into a single `MultiFab` per AMReX level.
     * Vector fields with three components are written with `x`, `y`, and `z`
     * suffixes. If `include_ghosts` is true, the output boxes are grown by the
     * source ghost-cell width and a matching output geometry is constructed.
     *
     * @param step Integer output step.
     * @param time Simulation time stored in the plotfile header.
     * @param include_ghosts If true, include ghost cells in the output.
     */
    void writeAMReXPlotfile(int step,
                            amrex::Real time = 0.0,
                            bool include_ghosts = false)
    {
        using namespace amrex;
        namespace fs = std::filesystem;

        if (fields.empty()) {
            amrex::Print() << "No fields to write.\n";
            return;
        }

        /*
         * Determine how many AMReX levels are represented by the registered
         * fields.
         */
        int finest_level = 0;

        for (auto& f : fields) {
            finest_level = std::max(finest_level, f.level);
        }

        int nlevels = finest_level + 1;

        /*
         * Flatten all registered component ranges into one global variable list.
         */
        Vector<std::string> var_names;
        int ncomp_total = 0;

        for (auto& f : fields) {
            if (!f.mf) continue;

            int ncomp_src = f.i_max - f.i_min;

            for (int c = f.i_min; c < f.i_max; ++c) {
                std::ostringstream ss;

                ss << f.name;

                if (f.i_max - f.i_min > 1) {
                    if (c - f.i_min == 0) ss << 'x';
                    if (c - f.i_min == 1) ss << 'y';
                    if (c - f.i_min == 2) ss << 'z';
                }

                var_names.push_back(ss.str());
            }

            ncomp_total += ncomp_src;
        }

        /*
         * Allocate the combined output MultiFab for each level. The first field
         * registered at a level defines that level's layout and geometry.
         */
        Vector<std::unique_ptr<MultiFab>> combined(nlevels);
        Vector<Geometry> geoms(nlevels);

        for (int lev = 0; lev < nlevels; ++lev) {
            auto it = std::find_if(fields.begin(), fields.end(),
                [lev](const auto& f) { return f.level == lev; });

            if (it == fields.end()) {
                amrex::Print() << "Missing data for level " << lev << "\n";
                continue;
            }

            const MultiFab& ref_mf = *(it->mf);
            const BoxArray& ba_in = ref_mf.boxArray();
            const DistributionMapping& dm = ref_mf.DistributionMap();
            const Geometry& geom_in = *(it->geom);
            const Real gsc = it->grid_scale;

            BoxArray ba_new = ba_in;

            if (include_ghosts) {
                int ng = ref_mf.nGrow();

                ba_new.grow(ng);
                geoms[lev] = MakeGeometryFromBoxArray(geom_in, ba_in, ba_new);
            } else {
                geoms[lev] = geom_in;
            }

            /*
             * Apply coordinate scaling to the output geometry only. Field data
             * scaling is handled separately when source components are copied.
             */
            if (apply_unit_scaling && gsc != 1.) {
                const auto plo = geoms[lev].ProbLoArray();
                const auto phi = geoms[lev].ProbHiArray();

                amrex::RealBox rb(
                    { plo[0] * gsc, plo[1] * gsc, plo[2] * gsc },
                    { phi[0] * gsc, phi[1] * gsc, phi[2] * gsc }
                );

                std::array<int, AMREX_SPACEDIM> is_per = {
                    geoms[lev].isPeriodic(0),
                    geoms[lev].isPeriodic(1),
                    geoms[lev].isPeriodic(2)
                };

                geoms[lev] = amrex::Geometry(geoms[lev].Domain(),
                                             &rb,
                                             amrex::CoordSys::cartesian,
                                             is_per.data());
            }

            combined[lev] = std::make_unique<amrex::MultiFab>(
                ba_new, dm, ncomp_total, 0
            );
        }

        /*
         * Copy each registered field into its assigned component range in the
         * combined output field.
         */
        for (int lev = 0; lev < nlevels; ++lev) {
            MultiFab& dst = *combined[lev];
            int comp_offset = 0;

            for (auto& f : fields) {
                if (!f.mf) continue;
                if (f.level != lev) continue;

                const MultiFab& full_src = *(f.mf);
                Real scale = f.scale;

                MultiFab src(full_src.boxArray(),
                             full_src.DistributionMap(),
                             f.i_max - f.i_min,
                             full_src.nGrow());

                src.setVal(0);

                MultiFab::Copy(src,
                               full_src,
                               f.i_min,
                               0,
                               f.i_max - f.i_min,
                               full_src.nGrow());

                if (apply_unit_scaling && scale != 1.) {
                    src.mult(scale, src.nGrow());
                }

                int ncomp_src = src.nComp();

                if (include_ghosts) {
                    dst.ParallelCopy(src,
                                     0,
                                     comp_offset,
                                     ncomp_src,
                                     src.nGrowVect(),
                                     dst.nGrowVect(),
                                     geoms[lev].periodicity());
                } else {
                    MultiFab::Copy(dst, src, 0, comp_offset, ncomp_src, 0);
                }

                comp_offset += ncomp_src;
            }
        }

        /*
         * AMReX expects one MultiFab pointer per level.
         */
        Vector<const MultiFab*> mf_ptrs(nlevels, nullptr);

        for (int lev = 0; lev < nlevels; ++lev) {
            mf_ptrs[lev] = combined[lev].get();
        }

        std::ostringstream folder;

        folder << output_dir
               << "/outputs"
               << database_name
               << "/plt"
               << std::setw(5)
               << std::setfill('0')
               << step;

        /*
         * The current solver does not use AMR refinement in these outputs, but
         * AMReX's multilevel writer still requires level step and refinement
         * ratio arrays.
         */
        Vector<int> level_steps(nlevels, step);
        Vector<IntVect> ref_ratios(std::max(0, nlevels - 1), IntVect(2));

        amrex::WriteMultiLevelPlotfile(folder.str(),
                                       nlevels,
                                       mf_ptrs,
                                       var_names,
                                       geoms,
                                       time,
                                       level_steps,
                                       ref_ratios);

#if PRINT_WHEN_SAVING_DATA
        amrex::Print() << "AMReX plotfile written with "
                       << ncomp_total
                       << " components and "
                       << nlevels
                       << " levels: "
                       << folder.str()
                       << (include_ghosts ? " (with ghosts)" : "")
                       << "\n";
#endif
    }

    /**
     * @brief Generate `.visit` time-series files for AMReX plotfile outputs.
     *
     * This routine scans the `outputs/` directory, finds AMReX plotfile
     * `Header` files for each output series, and writes one `.visit` file per
     * series. Paths are written relative to the directory containing the
     * `.visit` file.
     */
    void WriteAMReXMasterHeaderList() const
    {
        namespace fs = std::filesystem;

        const fs::path prefix = fs::path(output_dir);

        /*
         * Accept either a path that already ends in outputs/ or a path whose
         * child directory is outputs/.
         */
        fs::path outputs_dir;

        if (prefix.filename() == "outputs") {
            outputs_dir = prefix;
        } else if (fs::exists(prefix / "outputs") &&
                   fs::is_directory(prefix / "outputs")) {
            outputs_dir = prefix / "outputs";
        } else {
            std::cerr << "Error: could not locate an 'outputs' directory from m_plot_prefix = "
                      << prefix << "\n";
            return;
        }

        const fs::path visit_dir = outputs_dir.parent_path();

        if (!fs::exists(outputs_dir) || !fs::is_directory(outputs_dir)) {
            std::cerr << "Error: outputs directory does not exist: "
                      << outputs_dir << "\n";
            return;
        }

        bool wrote_any = false;

        /*
         * Each immediate child of outputs/ is treated as one plotfile series,
         * for example outputs/state or outputs/pml.
         */
        for (const auto& top_entry : fs::directory_iterator(outputs_dir))
        {
            if (!top_entry.is_directory()) continue;

            const fs::path series_dir = top_entry.path();
            const std::string series_name = series_dir.filename().string();

            std::vector<fs::path> headers;

            /*
             * Recursively collect AMReX plotfile Header files for this series.
             */
            for (const auto& entry : fs::recursive_directory_iterator(series_dir))
            {
                if (!entry.is_regular_file()) continue;
                if (entry.path().filename() != "Header") continue;

                const fs::path parent_dir = entry.path().parent_path();
                const std::string parent_name = parent_dir.filename().string();

                if (parent_name.rfind("plt", 0) != 0) continue;

                bool in_old_dir = false;

                for (const auto& part : entry.path())
                {
                    if (part.string().find("old") != std::string::npos) {
                        in_old_dir = true;
                        break;
                    }
                }

                if (in_old_dir) continue;

                headers.push_back(entry.path());
            }

            if (headers.empty()) {
                std::cout << "No plotfile headers found under "
                          << series_dir << "\n";
                continue;
            }

            std::sort(headers.begin(), headers.end());

            const fs::path visit_file = visit_dir / (series_name + ".visit");
            std::ofstream ofs(visit_file);

            if (!ofs.is_open()) {
                std::cerr << "Error: could not open "
                          << visit_file
                          << " for writing.\n";
                continue;
            }

            for (const auto& h : headers)
            {
                fs::path rel = fs::relative(h, visit_dir);
                ofs << rel.generic_string() << "\n";
            }

            ofs.close();

#if PRINT_WHEN_SAVING_DATA
            std::cout << "Wrote "
                      << visit_file
                      << " with "
                      << headers.size()
                      << " entries from "
                      << series_dir
                      << "\n";
#endif
            wrote_any = true;
        }

        if (!wrote_any) {
            std::cerr << "No AMReX plotfile Header files found under "
                      << outputs_dir << "\n";
        }
    }

    /**
     * @brief Coordinate-component indices used by some output helpers.
     */
    enum { X = 0, Y = 1, Z = 2 };

private:
    /**
     * @brief Output database/series suffix.
     */
    std::string database_name = "";

    /**
     * @brief Registered fields written by this manager.
     */
    std::vector<FieldData> fields;

    /**
     * @brief True when field and grid output scaling should be applied.
     */
    bool apply_unit_scaling = false;
};

#endif