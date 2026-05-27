#ifdef USE_IMAS
#include "io/ImasReader.H"
#include "io/Config.H"
#include "IMASH.H"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <string>

namespace
{
    inline int idx2(int iR, int iZ, int nR)
    {
        /*
         * IMAS-derived axisymmetric arrays are stored as a flattened
         * two-dimensional R-Z grid with R as the fastest index.
         */
        return iR + nR * iZ;
    }

    double pick_reference_frequency(const Config& cfg)
    {
        /*
         * The IMAS-driven computational grid is sized from a target number of
         * points per vacuum wavelength. Prefer an explicitly supplied reference
         * frequency; otherwise use the largest enabled source frequency.
         */
        if (cfg.reference_frequency > 0.0) {
            return cfg.reference_frequency;
        }

        double fmax = -1.0;

        for (const auto& s : cfg.sources) {
            if (s.enabled && s.frequency > fmax) {
                fmax = s.frequency;
            }
        }

        if (fmax > 0.0) {
            return fmax;
        }

        throw std::runtime_error(
            "BuildImashSetup: reference_frequency must be > 0, "
            "or at least one enabled source must have frequency > 0."
        );
    }

    int cells_from_extent(double L, double d_target)
    {
        /*
         * Convert a physical extent into an integer number of cells while
         * guaranteeing at least one cell. This is used for reduced-dimensional
         * cases as well as full 3-D domains.
         */
        if (L <= 0.0 || d_target <= 0.0) {
            return 1;
        }

        return std::max(1, static_cast<int>(std::ceil(L / d_target)));
    }

    double compute_dt_max_cartesian(double dx, double dy, double dz)
    {
        /*
         * Cartesian electromagnetic CFL limit for a three-dimensional FDTD
         * stencil. Reduced dimensions are represented by one-cell periodic
         * directions, so their cell size still appears in this estimate.
         */
        return 1.0 / (c_0 * std::sqrt(1.0 / (dx * dx)
                                    + 1.0 / (dy * dy)
                                    + 1.0 / (dz * dz)));
    }

    int find_bracket(const std::vector<double>& x, double xq)
    {
        /*
         * Return the lower index of the interval containing xq. A return value
         * of -1 indicates that the query lies outside the tabulated grid.
         */
        if (x.empty() || xq < x.front() || xq > x.back()) {
            return -1;
        }

        auto it = std::upper_bound(x.begin(), x.end(), xq);

        if (it == x.begin()) {
            return 0;
        }

        if (it == x.end()) {
            return static_cast<int>(x.size()) - 2;
        }

        return static_cast<int>((it - x.begin()) - 1);
    }

    double bilinear_on_grid(const AxisymmetricGrid2D& g,
                            const std::vector<double>& f,
                            double Rq,
                            double Zq,
                            bool& ok)
    {
        /*
         * Bilinear interpolation on the axisymmetric R-Z grid. The ok flag is
         * set only when the query lies inside a valid interpolation cell.
         */
        ok = false;

        if (g.nR < 2 || g.nZ < 2) {
            return 0.0;
        }

        const int iR = find_bracket(g.Rvals, Rq);
        const int iZ = find_bracket(g.Zvals, Zq);

        if (iR < 0 || iZ < 0 || iR >= g.nR - 1 || iZ >= g.nZ - 1) {
            return 0.0;
        }

        const double R0 = g.Rvals[iR];
        const double R1 = g.Rvals[iR + 1];
        const double Z0 = g.Zvals[iZ];
        const double Z1 = g.Zvals[iZ + 1];

        const double tx = (Rq - R0) / (R1 - R0);
        const double tz = (Zq - Z0) / (Z1 - Z0);

        const double f00 = f[idx2(iR,     iZ,     g.nR)];
        const double f10 = f[idx2(iR + 1, iZ,     g.nR)];
        const double f01 = f[idx2(iR,     iZ + 1, g.nR)];
        const double f11 = f[idx2(iR + 1, iZ + 1, g.nR)];

        ok = true;

        return (1.0 - tx) * (1.0 - tz) * f00
             + tx         * (1.0 - tz) * f10
             + (1.0 - tx) * tz         * f01
             + tx         * tz         * f11;
    }

    inline int clamp_int(int v, int lo, int hi)
    {
        return std::max(lo, std::min(v, hi));
    }

    double cubic_lagrange_1d(double x,
                             double x0, double x1, double x2, double x3,
                             double f0, double f1, double f2, double f3)
    {
        /*
         * Four-point Lagrange interpolation in one dimension. This helper is
         * used twice by bicubic_on_grid(): first in R, then in Z.
         */
        const double L0 =
            ((x - x1) * (x - x2) * (x - x3)) /
            ((x0 - x1) * (x0 - x2) * (x0 - x3));

        const double L1 =
            ((x - x0) * (x - x2) * (x - x3)) /
            ((x1 - x0) * (x1 - x2) * (x1 - x3));

        const double L2 =
            ((x - x0) * (x - x1) * (x - x3)) /
            ((x2 - x0) * (x2 - x1) * (x2 - x3));

        const double L3 =
            ((x - x0) * (x - x1) * (x - x2)) /
            ((x3 - x0) * (x3 - x1) * (x3 - x2));

        return f0 * L0 + f1 * L1 + f2 * L2 + f3 * L3;
    }

    double bicubic_on_grid(const AxisymmetricGrid2D& g,
                           const std::vector<double>& f,
                           double Rq,
                           double Zq,
                           bool& ok)
    {
        /*
         * Tensor-product cubic interpolation on the R-Z grid. Near boundaries,
         * or on grids too small to support a four-point stencil in both
         * directions, the routine falls back to bilinear interpolation.
         */
        ok = false;

        if (g.nR < 4 || g.nZ < 4) {
            return bilinear_on_grid(g, f, Rq, Zq, ok);
        }

        const int iR = find_bracket(g.Rvals, Rq);
        const int iZ = find_bracket(g.Zvals, Zq);

        if (iR < 0 || iZ < 0 || iR >= g.nR - 1 || iZ >= g.nZ - 1) {
            return 0.0;
        }

        const int iR0 = iR - 1;
        const int iZ0 = iZ - 1;

        if (iR0 < 0 || iR0 + 3 >= g.nR ||
            iZ0 < 0 || iZ0 + 3 >= g.nZ) {
            return bilinear_on_grid(g, f, Rq, Zq, ok);
        }

        const int r0 = iR0;
        const int r1 = iR0 + 1;
        const int r2 = iR0 + 2;
        const int r3 = iR0 + 3;

        const int z0 = iZ0;
        const int z1 = iZ0 + 1;
        const int z2 = iZ0 + 2;
        const int z3 = iZ0 + 3;

        const double R0 = g.Rvals[r0];
        const double R1 = g.Rvals[r1];
        const double R2 = g.Rvals[r2];
        const double R3 = g.Rvals[r3];

        const double Z0 = g.Zvals[z0];
        const double Z1 = g.Zvals[z1];
        const double Z2 = g.Zvals[z2];
        const double Z3 = g.Zvals[z3];

        const double q0 = cubic_lagrange_1d(
            Rq,
            R0, R1, R2, R3,
            f[idx2(r0, z0, g.nR)],
            f[idx2(r1, z0, g.nR)],
            f[idx2(r2, z0, g.nR)],
            f[idx2(r3, z0, g.nR)]
        );

        const double q1 = cubic_lagrange_1d(
            Rq,
            R0, R1, R2, R3,
            f[idx2(r0, z1, g.nR)],
            f[idx2(r1, z1, g.nR)],
            f[idx2(r2, z1, g.nR)],
            f[idx2(r3, z1, g.nR)]
        );

        const double q2 = cubic_lagrange_1d(
            Rq,
            R0, R1, R2, R3,
            f[idx2(r0, z2, g.nR)],
            f[idx2(r1, z2, g.nR)],
            f[idx2(r2, z2, g.nR)],
            f[idx2(r3, z2, g.nR)]
        );

        const double q3 = cubic_lagrange_1d(
            Rq,
            R0, R1, R2, R3,
            f[idx2(r0, z3, g.nR)],
            f[idx2(r1, z3, g.nR)],
            f[idx2(r2, z3, g.nR)],
            f[idx2(r3, z3, g.nR)]
        );

        ok = true;

        return cubic_lagrange_1d(Zq,
                                 Z0, Z1, Z2, Z3,
                                 q0, q1, q2, q3);
    }

    double interpolation_on_grid(const AxisymmetricGrid2D& g,
                                 const std::vector<double>& f,
                                 double Rq,
                                 double Zq,
                                 bool& ok)
    {
        /*
         * Compile-time switch between linear and higher-order interpolation.
         * The macro value is expected to be provided by the build/configuration
         * system.
         */
    #if IMAS_DATA_INTERPOLATION == linear
        return bilinear_on_grid(g, f, Rq, Zq, ok);
    #else
        return bicubic_on_grid(g, f, Rq, Zq, ok);
    #endif
    }
}

ImashSetup BuildImashSetup(const Config& cfg)
{
    /*
     * Build the complete IMAS-backed setup used by the FDTD run. This routine
     * loads equilibrium/plasma data, determines the computational box, chooses
     * grid resolution from points per wavelength, and computes a CFL-limited
     * timestep.
     */
    ImashSetup out;

    if (!cfg.use_imas) {
        throw std::runtime_error("BuildImashSetup: cfg.use_imas must be enabled.");
    }

    if (cfg.imas_path.empty()) {
        throw std::runtime_error("BuildImashSetup: imas_path must be provided.");
    }

    if (cfg.points_per_wavelength <= 0.0) {
        throw std::runtime_error("BuildImashSetup: points_per_wavelength must be positive.");
    }

    const std::string uri = "imas:hdf5?path=" + cfg.imas_path;

    IdsNs::IDS ids;
    ids.open(uri, OPEN_PULSE);

    /*
     * Extract the axisymmetric plasma profiles and optional magnetic data using
     * the requested IMAS time and merge/smoothing options.
     */
    out.plasma = ExtractPlasmaData(ids,
                                   cfg.imas_time,
                                   cfg.load_core != 0,
                                   cfg.load_bfield != 0,
                                   cfg.load_edge_on_grid != 0,
                                   cfg.merge_density != 0,
                                   cfg.merge_temperature != 0,
                                   cfg.smooth_rho,
                                   cfg.smooth_theta);

    if (out.plasma.empty()) {
        throw std::runtime_error("BuildImashSetup: ExtractPlasmaData returned an empty PlasmaData.");
    }

    out.axis = ExtractMagneticAxis(ids, cfg.imas_time);

    if (cfg.load_wall) {
        out.wall = ExtractWall2D(ids, cfg.imas_time);
    } else {
        out.wall.clear();
    }

    /*
     * Store metadata that is useful for diagnostics, logging, and later setup
     * stages without requiring the original IDS object.
     */
    out.plasma_nR = out.plasma.grid.nR;
    out.plasma_nZ = out.plasma.grid.nZ;
    out.magnetic_axis_R = out.axis.R;
    out.magnetic_axis_Z = out.axis.Z;
    out.n_wall_contours = static_cast<int>(out.wall.rz_contours.size());
    out.n_wall_meshes   = static_cast<int>(out.wall.surface_meshes.size());

    /*
     * Use the reference electromagnetic wavelength to determine the target grid
     * spacing.
     */
    const double freq = pick_reference_frequency(cfg);
    const double lambda0 = c_0 / freq;
    const double d_target = lambda0 / cfg.points_per_wavelength;

    const double Rmin_eq = out.plasma.grid.Rvals.front();
    const double Rmax_eq = out.plasma.grid.Rvals.back();
    const double Zmin_eq = out.plasma.grid.Zvals.front();
    const double Zmax_eq = out.plasma.grid.Zvals.back();

    const double Rmin_box = Rmin_eq - cfg.pad_R;
    const double Rmax_box = Rmax_eq + cfg.pad_R;
    const double Zmin_box = Zmin_eq - cfg.pad_Z;
    const double Zmax_box = Zmax_eq + cfg.pad_Z;

    /*
     * The AMReX x coordinate represents major radius R. The AMReX z coordinate
     * represents vertical coordinate Z. The y coordinate is either a one-cell
     * periodic direction for 2-D runs or a finite toroidal wedge for 3-D runs.
     */
    out.prob_lo[0] = Rmin_box;
    out.prob_hi[0] = Rmax_box;
    out.prob_lo[2] = Zmin_box;
    out.prob_hi[2] = Zmax_box;

    const bool use_user_R =
        std::abs(cfg.R_max - cfg.R_min) > TINY;

    const bool use_user_Z =
        std::abs(cfg.Z_max - cfg.Z_min) > TINY;

    const bool collapsed_user_Z =
        std::abs(cfg.Z_max - cfg.Z_min) <= TINY;

    /*
     * User-provided R bounds override the IMAS extent if they define a nonzero
     * interval.
     */
    if (use_user_R) {
        out.prob_lo[0] = std::min(cfg.R_min, cfg.R_max);
        out.prob_hi[0] = std::max(cfg.R_min, cfg.R_max);
    }

    /*
     * A nonzero user Z interval defines the physical Z range. A collapsed Z
     * interval is interpreted as a reduced-dimensional run with one periodic
     * cell in z centered on the requested Z position.
     */
    if (use_user_Z) {
        out.prob_lo[2] = std::min(cfg.Z_min, cfg.Z_max);
        out.prob_hi[2] = std::max(cfg.Z_min, cfg.Z_max);
    } else if (collapsed_user_Z) {
        out.prob_lo[2] = cfg.Z_min - d_target / 2.0;
        out.prob_hi[2] = cfg.Z_min + d_target / 2.0;
    }

    /*
     * Clamp R to the available IMAS domain including padding. For collapsed-Z
     * runs, do not clamp the artificial one-cell z thickness because it is a
     * computational thickness rather than an IMAS sampling interval.
     */
    if (out.prob_lo[0] < Rmin_box) out.prob_lo[0] = Rmin_box;
    if (out.prob_hi[0] > Rmax_box) out.prob_hi[0] = Rmax_box;

    if (!collapsed_user_Z) {
        if (out.prob_lo[2] < Zmin_box) out.prob_lo[2] = Zmin_box;
        if (out.prob_hi[2] > Zmax_box) out.prob_hi[2] = Zmax_box;
    }

    out.n_cell[0] = cells_from_extent(out.prob_hi[0] - out.prob_lo[0],
                                      d_target);

    if (collapsed_user_Z) {
        out.n_cell[2] = 1;
    } else {
        out.n_cell[2] = cells_from_extent(out.prob_hi[2] - out.prob_lo[2],
                                          d_target);
    }

    if (cfg.toroidal_phi_deg <= 0.0) {
        /*
         * Axisymmetric or R-Z run: use a single periodic y cell to keep the
         * Cartesian FDTD data structures three-dimensional while suppressing
         * variation in the toroidal direction.
         */
        out.is_3d = 0;

        out.prob_lo[1] = -0.5 * d_target;
        out.prob_hi[1] =  0.5 * d_target;
        out.n_cell[1] = 1;

        out.is_periodic[0] = 0;
        out.is_periodic[1] = 1;
        out.is_periodic[2] = 0;
    } else {
        /*
         * 3-D toroidal wedge. The Cartesian y extent is chosen from the outer
         * major radius and the requested toroidal angular width.
         */
        out.is_3d = 1;

        const double phi_rad = cfg.toroidal_phi_deg * M_PI / 180.0;
        const double y_half  = Rmax_box * std::sin(0.5 * phi_rad);

        out.prob_lo[1] = -y_half;
        out.prob_hi[1] =  y_half;
        out.n_cell[1] = cells_from_extent(out.prob_hi[1] - out.prob_lo[1],
                                          d_target);

        out.is_periodic[0] = 0;
        out.is_periodic[1] = 0;
        out.is_periodic[2] = 0;
    }

    /*
     * A collapsed Z interval represents a one-cell periodic z direction.
     */
    if (collapsed_user_Z) {
        out.is_periodic[2] = 1;
    }

    out.dx = (out.prob_hi[0] - out.prob_lo[0])
           / static_cast<double>(out.n_cell[0]);

    out.dy = (out.prob_hi[1] - out.prob_lo[1])
           / static_cast<double>(out.n_cell[1]);

    out.dz = (out.prob_hi[2] - out.prob_lo[2])
           / static_cast<double>(out.n_cell[2]);

    out.dt_max = compute_dt_max_cartesian(out.dx, out.dy, out.dz);
    out.dt     = cfg.cfl * out.dt_max;

    return out;
}

ImashSample SampleImashCartesian(const ImashSetup& setup,
                                 double x,
                                 double y,
                                 double z)
{
    /*
     * Sample axisymmetric IMAS data at a Cartesian point. The Cartesian point is
     * first converted to cylindrical coordinates, the scalar and cylindrical
     * magnetic components are interpolated on the R-Z grid, and the magnetic
     * field is then rotated back to Cartesian components.
     */
    ImashSample s;

    s.R   = std::sqrt(x * x + y * y);
    s.phi = std::atan2(y, x);
    s.Z   = z;

    const AxisymmetricGrid2D& g = setup.plasma.grid;

    const std::vector<double>& ne_h   = setup.plasma.ne.values;
    const std::vector<double>& Te_h   = setup.plasma.Te.values;
    const std::vector<double>& br_h   = setup.plasma.B.BR;
    const std::vector<double>& bz_h   = setup.plasma.B.BZ;
    const std::vector<double>& bphi_h = setup.plasma.B.Bphi;

    bool ok_ne   = false;
    bool ok_Te   = false;
    bool ok_br   = false;
    bool ok_bz   = false;
    bool ok_bphi = false;

    s.ne = interpolation_on_grid(g, ne_h, s.R, s.Z, ok_ne);
    s.Te = interpolation_on_grid(g, Te_h, s.R, s.Z, ok_Te);

    const double BR   = interpolation_on_grid(g, br_h,   s.R, s.Z, ok_br);
    const double BZ   = interpolation_on_grid(g, bz_h,   s.R, s.Z, ok_bz);
    const double Bphi = interpolation_on_grid(g, bphi_h, s.R, s.Z, ok_bphi);

    /*
     * The sample is considered valid only if every required plasma and magnetic
     * quantity was interpolated from inside the available R-Z grid.
     */
    s.inside_grid = ok_ne && ok_Te && ok_br && ok_bz && ok_bphi;

    const double cphi = std::cos(s.phi);
    const double sphi = std::sin(s.phi);

    /*
     * Cylindrical-to-Cartesian magnetic-field transform:
     *
     *   e_R   = (cos(phi), sin(phi), 0)
     *   e_phi = (-sin(phi), cos(phi), 0)
     *   e_Z   = (0, 0, 1)
     */
    s.Bx = BR * cphi - Bphi * sphi;
    s.By = BR * sphi + Bphi * cphi;
    s.Bz = BZ;

    return s;
}
#endif