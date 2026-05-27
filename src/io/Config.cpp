#include "io/Config.H"

#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <vector>

#include "utils/Constants.H"

void SourceConfig::remove_dimension()
{
    /*
     * Convert source parameters from physical units into the normalized units
     * used internally by the solver. The source shape factors are intentionally
     * not scaled here because they are interpreted as dimensionless stretch
     * factors rather than physical lengths.
     */
    position[0] /= L_0;
    position[1] /= L_0;
    position[2] /= L_0;

    width /= L_0;
    frequency /= nu_0;

    t_on   /= t_0;
    t_off  /= t_0;
    t_rise /= t_0;

    /*
     * The amplitude normalization depends on which field the source deposits
     * into.
     */
    if (type == "J" || type == "j") {
        amplitude /= J_0;
    } else if (type == "E" || type == "e") {
        amplitude /= E_0;
    } else if (type == "H" || type == "h") {
        amplitude /= H_0;
    }
}

void SourceConfig::add_dimension()
{
    /*
     * Restore source parameters from normalized solver units back to physical
     * units. This is the inverse of remove_dimension().
     */
    position[0] *= L_0;
    position[1] *= L_0;
    position[2] *= L_0;

    width *= L_0;
    frequency *= nu_0;

    t_on   *= t_0;
    t_off  *= t_0;
    t_rise *= t_0;

    if (type == "J" || type == "j") {
        amplitude *= J_0;
    } else if (type == "E" || type == "e") {
        amplitude *= E_0;
    } else if (type == "H" || type == "h") {
        amplitude *= H_0;
    }
}

void Config::read()
{
    /*
     * Read global initialization options. AMReX ParmParse leaves existing
     * defaults unchanged when a key is absent, so the Config defaults remain the
     * source of truth for unspecified inputs.
     */
    {
        amrex::ParmParse pp("init");

        pp.query("max_grid_size", max_grid_size);
        pp.query("nghost", nghost);
        pp.query("pml_thickness", pml_thickness);
        pp.query("cfl", cfl);
        pp.query("final_time", final_time);
        pp.query("max_step", max_step);
        pp.query("number_of_outputs", number_of_outputs);
        pp.query("use_pml", use_pml);
        pp.query("pml_sigma", pml_sigma);
        pp.query("plasma_sigma", plasma_sigma);
        pp.query("use_plasma", use_plasma);
        pp.query("use_e_cleaner", use_e_cleaner);
        pp.query("use_h_cleaner", use_h_cleaner);
        pp.query("output_dir", output_dir);
        pp.query("nsources", nsources);
        pp.query("normalize_units", normalize_units);

        /*
         * The code treats number_of_outputs as the number of output intervals
         * plus the initial state.
         */
        number_of_outputs++;
    }

    /*
     * Read manually specified Cartesian domain data. These values are used when
     * IMAS setup is disabled.
     */
    {
        amrex::ParmParse pp("input");

        pp.query("scale", input_scale);

        std::vector<int> vi(3);
        std::vector<double> vr(3);

        if (pp.queryarr("n_cell", vi, 0, 3)) {
            for (int d = 0; d < 3; ++d) input_n_cell[d] = vi[d];
        }

        if (pp.queryarr("prob_lo", vr, 0, 3)) {
            for (int d = 0; d < 3; ++d) input_prob_lo[d] = vr[d];
        }

        if (pp.queryarr("prob_hi", vr, 0, 3)) {
            for (int d = 0; d < 3; ++d) input_prob_hi[d] = vr[d];
        }

        vi.assign(3, 0);

        if (pp.queryarr("is_periodic", vi, 0, 3)) {
            for (int d = 0; d < 3; ++d) input_is_periodic[d] = vi[d];
        }
    }

#ifdef USE_IMAS
    /*
     * Read IMAS-specific options only when the executable is built with IMAS
     * support. The validation phase later decides whether this branch is
     * actually active for the current input file.
     */
    {
        amrex::ParmParse pp("imas");

        pp.query("use_imas", use_imas);
        pp.query("imas_path", imas_path);
        pp.query("imas_time", imas_time);
        pp.query("load_core", load_core);
        pp.query("load_bfield", load_bfield);
        pp.query("load_edge_on_grid", load_edge_on_grid);
        pp.query("load_wall", load_wall);
        pp.query("merge_density", merge_density);
        pp.query("merge_temperature", merge_temperature);
        pp.query("smooth_rho", smooth_rho);
        pp.query("smooth_theta", smooth_theta);
        pp.query("toroidal_phi_deg", toroidal_phi_deg);
        pp.query("points_per_wavelength", points_per_wavelength);
        pp.query("pad_R", pad_R);
        pp.query("pad_Z", pad_Z);
        pp.query("reference_frequency", reference_frequency);
        pp.query("R_min", R_min);
        pp.query("R_max", R_max);
        pp.query("Z_min", Z_min);
        pp.query("Z_max", Z_max);
    }
#endif

    /*
     * Read source.N blocks. Each source is optional, but the vector is sized to
     * match init.nsources so validation can check consistency explicitly.
     */
    sources.clear();
    sources.resize(nsources);

    for (int i = 0; i < nsources; ++i)
    {
        SourceConfig& s = sources[i];

        amrex::ParmParse pp("source." + std::to_string(i + 1));

        pp.query("enabled", s.enabled);
        pp.query("type", s.type);

        std::vector<double> vr(3);

        if (pp.queryarr("position", vr, 0, 3)) {
            for (int d = 0; d < 3; ++d) s.position[d] = vr[d];
        }

        vr.assign(3, 0.0);

        if (pp.queryarr("direction", vr, 0, 3)) {
            for (int d = 0; d < 3; ++d) s.direction[d] = vr[d];
        }

        vr.assign(3, 1.0);

        if (pp.queryarr("shape", vr, 0, 3)) {
            for (int d = 0; d < 3; ++d) s.shape[d] = vr[d];
        }

        pp.query("amplitude", s.amplitude);
        pp.query("frequency", s.frequency);
        pp.query("phase", s.phase);
        pp.query("width", s.width);
        pp.query("profile", s.profile);
        pp.query("t_on", s.t_on);
        pp.query("t_off", s.t_off);
        pp.query("t_rise", s.t_rise);
        pp.query("supergaussian_order", s.supergaussian_order);
    }
}

void Config::validate() const
{
    /*
     * Validation is split into local checks, branch-specific checks, source
     * checks, and cross-checks between options. This keeps error messages
     * targeted and makes future validation rules easier to add.
     */
    validate_common();
    validate_branch();
    validate_sources();
    validate_cross_checks();
}

void Config::require(bool cond, const std::string& msg)
{
    if (!cond) {
        throw std::runtime_error("Config validation failed: " + msg);
    }
}

bool Config::is_binary_flag(int x)
{
    return x == 0 || x == 1;
}

double Config::norm2(const std::array<double, 3>& a)
{
    return a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
}

void Config::validate_common() const
{
    /*
     * These checks apply to both manually specified domains and IMAS-derived
     * domains.
     */
    require(nghost >= 0, "init.nghost must be >= 0");
    require(pml_thickness >= 0, "init.pml_thickness must be >= 0");

    require(cfl > 0.0, "init.cfl must be > 0");
    require(final_time >= 0.0, "init.final_time must be >= 0");
    require(max_step == -1 || max_step >= 0,
            "init.max_step must be -1 or >= 0");
    require(number_of_outputs > 0.0,
            "init.number_of_outputs must be > 0");

    require(is_binary_flag(use_pml), "init.use_pml must be 0 or 1");
    require(is_binary_flag(use_plasma), "init.use_plasma must be 0 or 1");
    require(is_binary_flag(use_e_cleaner),
            "init.use_e_cleaner must be 0 or 1");
    require(is_binary_flag(use_h_cleaner),
            "init.use_h_cleaner must be 0 or 1");

    require(!output_dir.empty(), "init.output_dir must not be empty");
    require(nsources >= 0, "init.nsources must be >= 0");
    require(static_cast<int>(sources.size()) == nsources,
            "sources.size() must match init.nsources");
}

void Config::validate_branch() const
{
    /*
     * Only one domain-definition branch is active: IMAS-derived setup or manual
     * input setup.
     */
    if (use_imas) {
        validate_imas();
    } else {
        validate_manual_input();
    }
}

void Config::validate_manual_input() const
{
    /*
     * Manual input requires an explicit Cartesian box, cell count, and
     * periodicity flag in each coordinate direction.
     */
    require(input_scale > 0.0,
            "input.scale must be > 0 when imas.use_imas = 0");

    for (int d = 0; d < 3; ++d) {
        {
            std::ostringstream oss;
            oss << "input.n_cell[" << d << "] must be > 0";
            require(input_n_cell[d] > 0, oss.str());
        }

        {
            std::ostringstream oss;
            oss << "input.prob_hi[" << d
                << "] must be > input.prob_lo[" << d << "]";
            require(input_prob_hi[d] > input_prob_lo[d], oss.str());
        }

        {
            std::ostringstream oss;
            oss << "input.is_periodic[" << d << "] must be 0 or 1";
            require(is_binary_flag(input_is_periodic[d]), oss.str());
        }
    }
}

void Config::validate_imas() const
{
#ifdef USE_IMAS
    /*
     * IMAS validation is compiled only into IMAS-enabled builds. In non-IMAS
     * builds this routine is intentionally empty.
     */
    require(!imas_path.empty(),
            "imas.imas_path must not be empty when imas.use_imas = 1");
    require(imas_time >= 0.0, "imas.imas_time must be >= 0");

    require(is_binary_flag(load_core), "imas.load_core must be 0 or 1");
    require(is_binary_flag(load_bfield), "imas.load_bfield must be 0 or 1");
    require(is_binary_flag(load_edge_on_grid),
            "imas.load_edge_on_grid must be 0 or 1");
    require(is_binary_flag(load_wall), "imas.load_wall must be 0 or 1");
    require(is_binary_flag(merge_density),
            "imas.merge_density must be 0 or 1");
    require(is_binary_flag(merge_temperature),
            "imas.merge_temperature must be 0 or 1");

    require(smooth_rho >= 0.0, "imas.smooth_rho must be >= 0");
    require(smooth_theta >= 0.0, "imas.smooth_theta must be >= 0");
    require(points_per_wavelength > 0.0,
            "imas.points_per_wavelength must be > 0");
    require(pad_R >= 0.0, "imas.pad_R must be >= 0");
    require(pad_Z >= 0.0, "imas.pad_Z must be >= 0");

    require(reference_frequency == -1.0 || reference_frequency > 0.0,
            "imas.reference_frequency must be -1 or > 0");

    require(R_max > R_min, "imas.R_max must be > imas.R_min");
    require(Z_max >= Z_min, "imas.Z_max must be >= imas.Z_min");

    require(toroidal_phi_deg >= 0.0,
            "imas.toroidal_phi_deg must be >= 0");
    require(toroidal_phi_deg <= 360.0,
            "imas.toroidal_phi_deg must be <= 360");
#endif
}

void Config::validate_sources() const
{
    for (int i = 0; i < nsources; ++i) {
        validate_one_source(sources[i], i + 1);
    }
}

void Config::validate_one_source(const SourceConfig& s, int idx) const
{
    /*
     * Disabled sources only need a valid enabled flag. All other source fields
     * are ignored for disabled entries.
     */
    {
        std::ostringstream oss;
        oss << "source." << idx << ".enabled must be 0 or 1";
        require(is_binary_flag(s.enabled), oss.str());
    }

    if (!s.enabled) return;

    {
        std::ostringstream oss;
        oss << "source." << idx << ".type must not be empty";
        require(!s.type.empty(), oss.str());
    }

    {
        std::ostringstream oss;
        oss << "source." << idx << ".frequency must be > 0";
        require(s.frequency > 0.0, oss.str());
    }

    {
        std::ostringstream oss;
        oss << "source." << idx << ".width must be > 0";
        require(s.width > 0.0, oss.str());
    }

    /*
     * Temporal envelope parameters are optional. If all three are provided, the
     * envelope must define nonnegative endpoints and a positive transition time.
     */
    const bool has_envelope =
        s.t_on != -1 && s.t_off != -1 && s.t_rise != -1;

    if (has_envelope) {
        {
            std::ostringstream oss;
            oss << "source." << idx << ".t_on must be >= 0";
            require(s.t_on >= 0.0, oss.str());
        }

        {
            std::ostringstream oss;
            oss << "source." << idx << ".t_off must be >= 0";
            require(s.t_off >= 0.0, oss.str());
        }

        {
            std::ostringstream oss;
            oss << "source." << idx << ".t_rise must be > 0";
            require(s.t_rise > 0.0, oss.str());
        }
    }

    {
        std::ostringstream oss;
        oss << "source." << idx << ".profile must not be empty";
        require(!s.profile.empty(), oss.str());
    }

    {
        std::ostringstream oss;
        oss << "source." << idx << ".supergaussian_order must be > 0";
        require(s.supergaussian_order > 0.0, oss.str());
    }

    {
        std::ostringstream oss;
        oss << "source." << idx << ".direction must be nonzero";
        require(norm2(s.direction) > 0.0, oss.str());
    }

    /*
     * For manual domains, source positions can be checked immediately against
     * the configured physical box. For IMAS domains, the final domain may be
     * constructed later from equilibrium data, so this check is skipped here.
     */
    if (!use_imas) {
        for (int d = 0; d < 3; ++d) {
            std::ostringstream oss;
            oss << "source." << idx
                << ".position[" << d
                << "] lies outside [input.prob_lo, input.prob_hi]";

            require(s.position[d] >= input_prob_lo[d]
                 && s.position[d] <= input_prob_hi[d],
                    oss.str());
        }
    }
}

int Config::count_enabled_sources() const
{
    int n = 0;

    for (const auto& s : sources) {
        if (s.enabled) ++n;
    }

    return n;
}

double Config::max_enabled_source_frequency() const
{
    double fmax = -1.0;

    for (const auto& s : sources) {
        if (s.enabled) {
            fmax = std::max(fmax, s.frequency);
        }
    }

    return fmax;
}

void Config::validate_cross_checks() const
{
    /*
     * Cross-checks compare options that are individually valid but incompatible
     * when used together.
     */
    validate_pml_vs_manual_domain();
    validate_imas_vs_sources();
}

void Config::validate_pml_vs_manual_domain() const
{
    /*
     * For manually specified domains, an explicitly cell-counted PML must fit
     * inside each non-periodic direction. Fractional PML thickness is handled by
     * the shell-generation routine and is not checked here.
     */
    if (!use_pml || use_imas) return;

    if (pml_thickness > 1) {
        for (int d = 0; d < 3; ++d) {
            if (input_is_periodic[d] == 0) {
                std::ostringstream oss;
                oss << "init.pml_thickness = " << pml_thickness
                    << " is too large for input.n_cell[" << d << "] = "
                    << input_n_cell[d]
                    << " (expect 2*pml_thickness < n_cell in non-periodic directions)";

                require(2 * pml_thickness < input_n_cell[d] * input_scale,
                        oss.str());
            }
        }
    }
}

void Config::validate_imas_vs_sources() const
{
#ifdef USE_IMAS
    /*
     * IMAS grid construction needs either an explicit reference frequency or at
     * least one enabled source frequency from which a wavelength can be inferred.
     */
    if (!use_imas) return;

    const int n_enabled = count_enabled_sources();
    const double fmax = max_enabled_source_frequency();

    require(n_enabled > 0,
            "at least one enabled source is required when imas.use_imas = 1");

    if (reference_frequency < 0.0) {
        require(fmax > 0.0,
                "imas.reference_frequency < 0 requires at least one enabled source with frequency > 0");
    }
#endif
}

void Config::remove_dimension()
{
    /*
     * Convert global configuration quantities from physical units into solver
     * normalized units. Sentinel values such as BIGGY-sized bounds are left
     * untouched.
     */
    final_time /= t_0;

    input_prob_lo[0] /= L_0;
    input_prob_lo[1] /= L_0;
    input_prob_lo[2] /= L_0;

    input_prob_hi[0] /= L_0;
    input_prob_hi[1] /= L_0;
    input_prob_hi[2] /= L_0;

    pad_R /= L_0;
    pad_Z /= L_0;

    if (reference_frequency > 0.0) {
        reference_frequency /= nu_0;
    }

    if (R_min > -BIGGY / 2.0 && R_min < BIGGY / 2.0) R_min /= L_0;
    if (R_max > -BIGGY / 2.0 && R_max < BIGGY / 2.0) R_max /= L_0;
    if (Z_min > -BIGGY / 2.0 && Z_min < BIGGY / 2.0) Z_min /= L_0;
    if (Z_max > -BIGGY / 2.0 && Z_max < BIGGY / 2.0) Z_max /= L_0;

    for (auto& s : sources) {
        s.remove_dimension();
    }
}

void Config::add_dimension()
{
    /*
     * Restore global configuration quantities from normalized solver units back
     * to physical units. Sentinel bounds are left untouched.
     */
    final_time *= t_0;

    input_prob_lo[0] *= L_0;
    input_prob_lo[1] *= L_0;
    input_prob_lo[2] *= L_0;

    input_prob_hi[0] *= L_0;
    input_prob_hi[1] *= L_0;
    input_prob_hi[2] *= L_0;

    pad_R *= L_0;
    pad_Z *= L_0;

    if (reference_frequency > 0.0) {
        reference_frequency *= nu_0;
    }

    if (R_min > -BIGGY / 2.0 && R_min < BIGGY / 2.0) R_min *= L_0;
    if (R_max > -BIGGY / 2.0 && R_max < BIGGY / 2.0) R_max *= L_0;
    if (Z_min > -BIGGY / 2.0 && Z_min < BIGGY / 2.0) Z_min *= L_0;
    if (Z_max > -BIGGY / 2.0 && Z_max < BIGGY / 2.0) Z_max *= L_0;

    for (auto& s : sources) {
        s.add_dimension();
    }
}