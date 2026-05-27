#ifndef POINT_SOURCE_H_
#define POINT_SOURCE_H_

#include "sources/Source.H"
#include "fdtd/EHJSolver.H"

#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

/**
 * @brief Localized electromagnetic source deposited into `J`, `E`, or `H`.
 *
 * `PointSource` injects a time-dependent source into one of the solver fields.
 * Despite the name, the source can have finite spatial extent through a
 * Gaussian, super-Gaussian, or Heaviside profile. The profile is centered at a
 * physical position and can be stretched independently along each coordinate
 * direction using the `shape` factors.
 *
 * The temporal waveform can be constant or sinusoidal. For oscillatory sources,
 * an optional smooth tanh turn-on/turn-off envelope can be applied to avoid
 * sharp startup and shutdown transients.
 *
 * The source stores a pointer to an existing electromagnetic solver. The solver
 * and geometry objects must outlive the source.
 */
class PointSource : public Source
{
public:
    using SourceType  = Source::SourceType;
    using ProfileType = Source::ProfileType;

    /**
     * @brief Construct a localized source attached to an electromagnetic solver.
     *
     * @param solver_ Solver whose fields receive the source contribution.
     * @param type_ Target field: current density `J`, electric field `E`, or
     * magnetic field `H`.
     * @param x_ Source-center x coordinate in physical/code geometry units.
     * @param y_ Source-center y coordinate in physical/code geometry units.
     * @param z_ Source-center z coordinate in physical/code geometry units.
     * @param amplitude_ Source amplitude before temporal and spatial weighting.
     * @param frequency_ Oscillation frequency. A nonpositive value gives a
     * non-oscillatory source.
     * @param phase_ Oscillation phase in degrees.
     * @param width_ Spatial profile width. For a Heaviside profile this is the
     * support radius.
     * @param profile_ Spatial profile used to distribute the source.
     * @param supergaussian_order_ Order used when `profile_` is
     * `ProfileType::SuperGaussian`.
     * @param direction_ Three-component source polarization/direction vector.
     * The vector is normalized internally.
     * @param geom_ Geometry descriptor used to convert physical coordinates to
     * cell indices.
     * @param shape_ Three-component anisotropic stretch factors for the source
     * profile.
     * @param t_on Start time for the smooth source envelope. Use `-1` together
     * with `t_off = -1` to disable the envelope.
     * @param t_off Stop time for the smooth source envelope. Use `-1` together
     * with `t_on = -1` to disable the envelope.
     * @param t_rise Characteristic rise/fall time of the tanh envelope.
     */
    PointSource(EHSolver& solver_,
                SourceType type_,
                amrex::Real x_, amrex::Real y_, amrex::Real z_,
                amrex::Real amplitude_,
                amrex::Real frequency_,
                amrex::Real phase_,
                amrex::Real width_,
                ProfileType profile_,
                amrex::Real supergaussian_order_,
                amrex::Vector<amrex::Real> const& direction_,
                const amrex::Geometry& geom_,
                amrex::Vector<amrex::Real> const& shape_,
                amrex::Real t_on = -1,
                amrex::Real t_off = -1,
                amrex::Real t_rise = -1);

    /**
     * @brief Evaluate and deposit the source at a given simulation time.
     *
     * The method computes the temporal source amplitude, applies the optional
     * turn-on/turn-off envelope, and deposits the result into the selected
     * field.
     *
     * @param time Current simulation time.
     * @param dt Timestep size. For direct `E` and `H` sources this is used to
     * scale the deposited field increment.
     */
    void update(amrex::Real time, amrex::Real dt = 1.0) override;

    /**
     * @brief Deposit the source into a specific vector field.
     *
     * The value is distributed over the source support using the configured
     * spatial profile and direction vector. The target `MultiFab` is assumed to
     * have at least three components.
     *
     * @param mf Target field receiving the source contribution.
     * @param val Scalar source amplitude after temporal weighting.
     */
    void setPoint(amrex::MultiFab& mf, amrex::Real val);

private:
    /**
     * @brief Solver whose fields receive the source contribution.
     *
     * The source does not own the solver.
     */
    EHSolver* EHsolver = nullptr;

    /**
     * @brief Field targeted by the source.
     */
    SourceType type;

    /**
     * @brief Geometry used to map source coordinates to grid cells.
     *
     * Stored by reference; the referenced object must outlive the source.
     */
    const amrex::Geometry& geom;

    /**
     * @brief Source-center x coordinate.
     */
    amrex::Real x0 = 0.0;

    /**
     * @brief Source-center y coordinate.
     */
    amrex::Real y0 = 0.0;

    /**
     * @brief Source-center z coordinate.
     */
    amrex::Real z0 = 0.0;

    /**
     * @brief Normalized three-component source direction.
     */
    amrex::GpuArray<amrex::Real, 3> dir;

    /**
     * @brief Cached lower corner of the computational domain.
     */
    amrex::GpuArray<amrex::Real, 3> prob_lo;

    /**
     * @brief Cached grid spacing in each coordinate direction.
     */
    amrex::GpuArray<amrex::Real, 3> dx;

    /**
     * @brief Anisotropic stretch factors for the spatial source profile.
     */
    amrex::GpuArray<amrex::Real, 3> shape;

    /**
     * @brief Source amplitude before temporal and spatial weighting.
     */
    amrex::Real amplitude = 0.0;

    /**
     * @brief Source oscillation frequency.
     */
    amrex::Real frequency = 0.0;

    /**
     * @brief Source oscillation phase in degrees.
     */
    amrex::Real phase = 0.0;

    /**
     * @brief Width of the spatial profile.
     */
    amrex::Real width = 1.0;

    /**
     * @brief Spatial profile used for deposition.
     */
    ProfileType profile = ProfileType::Gaussian;

    /**
     * @brief Super-Gaussian order used when `profile` is
     * `ProfileType::SuperGaussian`.
     */
    amrex::Real supergaussian_order = 2.0;

    /**
     * @brief Start time of the smooth temporal envelope.
     */
    amrex::Real t_start = -1;

    /**
     * @brief Stop time of the smooth temporal envelope.
     */
    amrex::Real t_stop = -1;

    /**
     * @brief Characteristic rise/fall time of the temporal envelope.
     */
    amrex::Real t_switch = -1;
};

#endif