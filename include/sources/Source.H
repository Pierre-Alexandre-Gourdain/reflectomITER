#ifndef SOURCE_H_
#define SOURCE_H_

#include <AMReX_REAL.H>

/**
 * @brief Abstract base class for time-dependent electromagnetic sources.
 *
 * `Source` defines the common interface used by objects that inject current,
 * electric-field, or magnetic-field contributions into an FDTD solver. Concrete
 * source implementations specify the spatial deposition rule and the temporal
 * waveform.
 *
 * Sources are updated by calling `update()` with the current simulation time
 * and a scaling parameter, typically the timestep or a timestep-dependent
 * factor.
 */
class Source {
public:
    /**
     * @brief Field quantity targeted by the source.
     *
     * A source can deposit into the current density `J`, the electric field
     * `E`, or the magnetic field `H`.
     */
    enum class SourceType {
        J,  ///< Deposit into the current-density field.
        E,  ///< Deposit directly into the electric field.
        H   ///< Deposit directly into the magnetic field.
    };

    /**
     * @brief Spatial profile used to distribute a localized source.
     *
     * Concrete source classes can use these profiles to control how the source
     * amplitude is distributed around its nominal center.
     */
    enum class ProfileType {
        Gaussian,       ///< Smooth Gaussian spatial profile.
        SuperGaussian,  ///< Higher-order Gaussian-like profile with flatter core.
        Heaviside       ///< Compact top-hat profile.
    };

    /**
     * @brief Destroy the source object.
     */
    virtual ~Source() = default;

    /**
     * @brief Apply the source contribution at a given simulation time.
     *
     * Concrete implementations use this method to evaluate the source waveform
     * and deposit it into the target field.
     *
     * @param time Current simulation time.
     * @param scale Source scaling parameter, commonly the timestep `dt`.
     */
    virtual void update(amrex::Real time, amrex::Real scale) = 0;

private:
    /**
     * @brief Target field selected for this source.
     *
     * This member is currently private and unused by the base interface. Derived
     * classes may store their own target-field state if they need direct access.
     */
    SourceType type;
};

#endif