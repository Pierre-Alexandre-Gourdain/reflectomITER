#ifdef USE_IMAS
#ifndef IMAS_READER_H_
#define IMAS_READER_H_

#include "io/Config.H"
#include "TokamakData.H"
#include "utils/Constants.H"

#include <array>
#include <vector>

/**
 * @brief Plasma and magnetic-field sample obtained from IMAS data at one
 * Cartesian point.
 *
 * `ImashSample` stores both the cylindrical coordinates used for sampling the
 * axisymmetric IMAS data and the Cartesian plasma/magnetic quantities returned
 * to the FDTD code.
 */
struct ImashSample
{
    /**
     * @brief Cylindrical major radius of the sample point.
     */
    double R = 0.0;

    /**
     * @brief Toroidal angle of the sample point.
     */
    double phi = 0.0;

    /**
     * @brief Vertical coordinate of the sample point.
     */
    double Z = 0.0;

    /**
     * @brief Interpolated electron density.
     */
    double ne = 0.0;

    /**
     * @brief Interpolated electron temperature.
     */
    double Te = 0.0;

    /**
     * @brief Cartesian x component of the interpolated magnetic field.
     */
    double Bx = 0.0;

    /**
     * @brief Cartesian y component of the interpolated magnetic field.
     */
    double By = 0.0;

    /**
     * @brief Cartesian z component of the interpolated magnetic field.
     */
    double Bz = 0.0;

    /**
     * @brief True if all requested quantities were sampled inside the available
     * IMAS interpolation grid.
     */
    bool inside_grid = false;
};

/**
 * @brief Complete IMAS-derived setup used to initialize an FDTD run.
 *
 * `ImashSetup` stores the extracted axisymmetric plasma data, magnetic-axis
 * metadata, optional wall data, and the Cartesian computational domain derived
 * from the IMAS data and the user configuration.
 *
 * The AMReX computational coordinates are interpreted as:
 *
 * - `x`: major-radius-like Cartesian coordinate,
 * - `y`: toroidal/cartesian transverse coordinate,
 * - `z`: vertical coordinate.
 *
 * For reduced-dimensional runs, one coordinate direction may be represented by
 * a single periodic cell.
 */
struct ImashSetup
{
    /**
     * @brief Extracted plasma profiles and magnetic-field data on an
     * axisymmetric R-Z grid.
     */
    PlasmaData plasma;

    /**
     * @brief Magnetic-axis location extracted from IMAS.
     */
    MagneticAxisPoint axis;

    /**
     * @brief Wall contours and surface meshes extracted from IMAS, when
     * requested.
     */
    WallData wall;

    /**
     * @brief Number of radial grid points in the extracted plasma grid.
     */
    int plasma_nR = 0;

    /**
     * @brief Number of vertical grid points in the extracted plasma grid.
     */
    int plasma_nZ = 0;

    /**
     * @brief Major-radius coordinate of the magnetic axis.
     */
    double magnetic_axis_R = 0.0;

    /**
     * @brief Vertical coordinate of the magnetic axis.
     */
    double magnetic_axis_Z = 0.0;

    /**
     * @brief Number of 2-D wall contours loaded from IMAS.
     */
    int n_wall_contours = 0;

    /**
     * @brief Number of wall surface meshes loaded from IMAS.
     */
    int n_wall_meshes = 0;

    /**
     * @brief Nonzero when the generated computational domain is a 3-D toroidal
     * wedge.
     */
    int is_3d = 0;

    /**
     * @brief Lower physical bounds of the generated Cartesian computational
     * domain.
     */
    std::array<double, 3> prob_lo{{0.0, 0.0, 0.0}};

    /**
     * @brief Upper physical bounds of the generated Cartesian computational
     * domain.
     */
    std::array<double, 3> prob_hi{{0.0, 0.0, 0.0}};

    /**
     * @brief Number of cells in each Cartesian coordinate direction.
     */
    std::array<int, 3> n_cell{{1, 1, 1}};

    /**
     * @brief Periodicity flags for the generated computational domain.
     */
    std::array<int, 3> is_periodic{{0, 0, 0}};

    /**
     * @brief Cell size in the x direction.
     */
    double dx = 0.0;

    /**
     * @brief Cell size in the y direction.
     */
    double dy = 0.0;

    /**
     * @brief Cell size in the z direction.
     */
    double dz = 0.0;

    /**
     * @brief Maximum stable timestep estimated from the Cartesian
     * electromagnetic CFL condition.
     */
    double dt_max = 0.0;

    /**
     * @brief Selected timestep after applying the configured CFL number.
     */
    double dt = 0.0;
};

/**
 * @brief Build an IMAS-derived computational setup from the runtime
 * configuration.
 *
 * This routine opens the configured IMAS data source, extracts plasma and
 * magnetic data, constructs the Cartesian computational domain, chooses the
 * number of cells from the requested points per wavelength, and computes the
 * corresponding CFL-limited timestep.
 *
 * @param cfg Runtime configuration.
 * @return Complete IMAS-derived setup for grid generation and plasma
 * initialization.
 */
ImashSetup BuildImashSetup(const Config& cfg);

/**
 * @brief Sample IMAS-derived plasma data at a Cartesian point.
 *
 * The point `(x, y, z)` is converted to cylindrical coordinates
 * `(R, phi, Z)`. Scalar plasma quantities and cylindrical magnetic-field
 * components are interpolated on the axisymmetric R-Z grid, and the magnetic
 * field is then transformed back to Cartesian components.
 *
 * @param setup IMAS-derived setup containing the plasma interpolation grid.
 * @param x Cartesian x coordinate.
 * @param y Cartesian y coordinate.
 * @param z Cartesian z coordinate.
 * @return Interpolated plasma and magnetic-field sample.
 */
ImashSample SampleImashCartesian(const ImashSetup& setup,
                                 double x,
                                 double y,
                                 double z);

#endif
#endif