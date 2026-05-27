#ifndef CONFIG_H_
#define CONFIG_H_

#include <array>
#include <string>
#include <vector>

#include "utils/Constants.H"

/**
 * @brief Runtime configuration for one electromagnetic source.
 *
 * `SourceConfig` stores the parameters read from a `source.N` input block. The
 * source can target the current density `J`, electric field `E`, or magnetic
 * field `H`, and can use a finite-width spatial profile.
 *
 * The unit-conversion methods convert source parameters between physical units
 * and the normalized units used internally by the solver.
 */
struct SourceConfig
{
    /**
     * @brief Nonzero if this source is active.
     */
    int enabled = 0;

    /**
     * @brief Target field type.
     *
     * Expected values are `"J"`, `"E"`, or `"H"`, case-insensitively.
     */
    std::string type = "J";

    /**
     * @brief Source-center position.
     */
    std::array<double, 3> position{{0.0, 0.0, 0.0}};

    /**
     * @brief Source polarization or deposition direction.
     */
    std::array<double, 3> direction{{0.0, 0.0, 0.0}};

    /**
     * @brief Dimensionless anisotropic stretch factors for the spatial profile.
     */
    std::array<double, 3> shape{{1.0, 1.0, 1.0}};

    /**
     * @brief Source amplitude before temporal and spatial weighting.
     */
    double amplitude = 0.0;

    /**
     * @brief Source oscillation frequency.
     */
    double frequency = 0.0;

    /**
     * @brief Source phase in degrees.
     */
    double phase = 0.0;

    /**
     * @brief Width of the spatial source profile.
     */
    double width = 0.0;

    /**
     * @brief Optional source turn-on time.
     *
     * A value of `-1` disables the smooth temporal envelope when used together
     * with `t_off = -1` and `t_rise = -1`.
     */
    double t_on = -1.0;

    /**
     * @brief Optional source turn-off time.
     */
    double t_off = -1.0;

    /**
     * @brief Optional source rise/fall time for the smooth tanh envelope.
     */
    double t_rise = -1.0;

    /**
     * @brief Spatial source profile name.
     *
     * Expected values are profile names such as `"gaussian"`,
     * `"supergaussian"`, or `"heaviside"`, depending on the source factory.
     */
    std::string profile = "gaussian";

    /**
     * @brief Order used by the super-Gaussian spatial profile.
     */
    double supergaussian_order = 4.0;

    /**
     * @brief Convert source parameters from physical units to normalized solver
     * units.
     */
    void remove_dimension();

    /**
     * @brief Restore source parameters from normalized solver units to physical
     * units.
     */
    void add_dimension();
};

/**
 * @brief Runtime configuration for the FDTD simulation.
 *
 * `Config` owns all user-facing input parameters read through AMReX
 * `ParmParse`. It contains global initialization options, manual-domain input,
 * optional IMAS setup options, and a list of source configurations.
 *
 * The configuration workflow is:
 *
 * 1. construct `Config`,
 * 2. call `read()` to populate values from the input file,
 * 3. call `validate()` to check consistency,
 * 4. optionally call `remove_dimension()` if normalized solver units are used.
 */
struct Config
{
    /**
     * @name init.* options
     * @{
     */

    /**
     * @brief Maximum AMReX box size.
     *
     * A negative value leaves the initial `BoxArray` unsplit.
     */
    int max_grid_size = -1;

    /**
     * @brief Number of ghost cells allocated for solver fields.
     */
    int nghost = 1;

    /**
     * @brief PML thickness.
     *
     * Values greater than one are interpreted as a number of cells. Values less
     * than or equal to one are interpreted as a fraction by the shell-generation
     * utility.
     */
    double pml_thickness = 0;

    /**
     * @brief CFL number used to compute the timestep.
     */
    double cfl = 0.9;

    /**
     * @brief Final simulation time.
     */
    double final_time = 0.0;

    /**
     * @brief Maximum number of timesteps.
     *
     * A value of `-1` disables the step-count limit.
     */
    int max_step = -1;

    /**
     * @brief Number of output intervals requested by the user.
     *
     * The implementation increments this after reading to account for the
     * initial output state.
     */
    int number_of_outputs = 10;

    /**
     * @brief Nonzero if the PML solver should be enabled.
     */
    int use_pml = 0;

    /**
     * @brief Nonzero if the ADE plasma model should be enabled.
     */
    int use_plasma = 0;

    /**
     * @brief Nonzero if electric divergence cleaning should be enabled.
     */
    int use_e_cleaner = 0;

    /**
     * @brief Nonzero if magnetic divergence cleaning should be enabled.
     */
    int use_h_cleaner = 0;

    /**
     * @brief Maximum PML damping strength.
     */
    double pml_sigma = 50.0;

    /**
     * @brief Conductivity scale used when coupling plasma damping to the PML
     * profile.
     */
    double plasma_sigma = 50.0;

    /**
     * @brief Root output directory.
     */
    std::string output_dir = ".";

    /**
     * @brief Number of source blocks to read.
     */
    int nsources = 0;

    /**
     * @brief Nonzero if the solver should normalize physical inputs.
     */
    int normalize_units = 0;

    /** @} */

    /**
     * @name input.* manual-domain options
     * @{
     */

    /**
     * @brief Integer refinement factor applied to manual cell counts.
     *
     * Directions with one cell are kept collapsed.
     */
    double input_scale = 1.0;

    /**
     * @brief Manual-domain cell count in each coordinate direction.
     */
    std::array<int, 3> input_n_cell{{0, 0, 0}};

    /**
     * @brief Lower physical bounds of the manual computational domain.
     */
    std::array<double, 3> input_prob_lo{{0.0, 0.0, 0.0}};

    /**
     * @brief Upper physical bounds of the manual computational domain.
     */
    std::array<double, 3> input_prob_hi{{0.0, 0.0, 0.0}};

    /**
     * @brief Periodicity flags for the manual computational domain.
     */
    std::array<int, 3> input_is_periodic{{0, 0, 0}};

    /** @} */

    /**
     * @name imas.* options
     * @{
     */

    /**
     * @brief Nonzero if the computational setup should be derived from IMAS.
     */
    int use_imas = 0;

    /**
     * @brief IMAS HDF5 path or data source path.
     */
    std::string imas_path = "";

    /**
     * @brief IMAS time slice to load.
     */
    double imas_time = 0.0;

    /**
     * @brief Nonzero if core plasma profiles should be loaded.
     */
    int load_core = 0;

    /**
     * @brief Nonzero if magnetic-field data should be loaded.
     */
    int load_bfield = 0;

    /**
     * @brief Nonzero if edge plasma data should be loaded on the equilibrium
     * grid.
     */
    int load_edge_on_grid = 0;

    /**
     * @brief Nonzero if wall data should be loaded.
     */
    int load_wall = 0;

    /**
     * @brief Nonzero if core and edge density data should be merged.
     */
    int merge_density = 0;

    /**
     * @brief Nonzero if core and edge temperature data should be merged.
     */
    int merge_temperature = 0;

    /**
     * @brief Radial smoothing parameter used when merging IMAS profiles.
     */
    double smooth_rho = 0.0;

    /**
     * @brief Angular smoothing parameter used when merging IMAS profiles.
     */
    double smooth_theta = 0.0;

    /**
     * @brief Toroidal wedge angle in degrees.
     *
     * A nonpositive value selects a reduced-dimensional axisymmetric setup.
     */
    double toroidal_phi_deg = 0.0;

    /**
     * @brief Target number of grid points per reference wavelength.
     */
    double points_per_wavelength = 15.0;

    /**
     * @brief Radial padding added around the IMAS plasma domain.
     */
    double pad_R = 0.0;

    /**
     * @brief Vertical padding added around the IMAS plasma domain.
     */
    double pad_Z = 0.0;

    /**
     * @brief Reference frequency used to define the grid spacing.
     *
     * A value of `-1` means the frequency should be inferred from enabled
     * sources.
     */
    double reference_frequency = -1.0;

    /**
     * @brief Optional lower radial bound for the IMAS-derived domain.
     */
    double R_min = 0.0;

    /**
     * @brief Optional upper radial bound for the IMAS-derived domain.
     */
    double R_max = BIGGY;

    /**
     * @brief Optional lower vertical bound for the IMAS-derived domain.
     */
    double Z_min = -BIGGY;

    /**
     * @brief Optional upper vertical bound for the IMAS-derived domain.
     */
    double Z_max = BIGGY;

    /** @} */

    /**
     * @brief Source configurations read from `source.N` input blocks.
     */
    std::vector<SourceConfig> sources;

    /**
     * @brief Read configuration values from AMReX `ParmParse`.
     *
     * Missing keys keep their default values.
     */
    void read();

    /**
     * @brief Validate the configuration for internal consistency.
     *
     * Throws `std::runtime_error` with a descriptive message if validation
     * fails.
     */
    void validate() const;

    /**
     * @brief Convert dimensional configuration values to normalized solver
     * units.
     */
    void remove_dimension();

    /**
     * @brief Restore configuration values from normalized solver units to
     * physical units.
     */
    void add_dimension();

private:
    /**
     * @brief Throw a validation error if a condition is false.
     *
     * @param cond Condition that must be true.
     * @param msg Error message appended to the validation failure prefix.
     */
    static void require(bool cond, const std::string& msg);

    /**
     * @brief Return true if an integer is a binary flag.
     *
     * @param x Value to test.
     * @return True if `x` is `0` or `1`.
     */
    static bool is_binary_flag(int x);

    /**
     * @brief Return the squared Euclidean norm of a three-component vector.
     *
     * @param a Vector to evaluate.
     * @return Squared norm.
     */
    static double norm2(const std::array<double, 3>& a);

    /**
     * @brief Validate options common to all setup branches.
     */
    void validate_common() const;

    /**
     * @brief Validate either the manual-input branch or the IMAS branch.
     */
    void validate_branch() const;

    /**
     * @brief Validate manually specified domain inputs.
     */
    void validate_manual_input() const;

    /**
     * @brief Validate IMAS-specific inputs.
     */
    void validate_imas() const;

    /**
     * @brief Validate all configured source blocks.
     */
    void validate_sources() const;

    /**
     * @brief Validate one source block.
     *
     * @param s Source configuration to validate.
     * @param idx One-based source index used in error messages.
     */
    void validate_one_source(const SourceConfig& s, int idx) const;

    /**
     * @brief Count enabled sources.
     *
     * @return Number of enabled source configurations.
     */
    int count_enabled_sources() const;

    /**
     * @brief Return the maximum frequency among enabled sources.
     *
     * @return Maximum enabled source frequency, or `-1` if none are enabled.
     */
    double max_enabled_source_frequency() const;

    /**
     * @brief Run validation checks involving multiple configuration sections.
     */
    void validate_cross_checks() const;

    /**
     * @brief Validate PML thickness against a manually specified domain.
     */
    void validate_pml_vs_manual_domain() const;

    /**
     * @brief Validate IMAS setup requirements involving enabled sources.
     */
    void validate_imas_vs_sources() const;
};

#endif