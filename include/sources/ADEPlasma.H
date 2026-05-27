#pragma once

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include "io/FieldOutputManager.H"

#ifdef USE_IMAS
struct ImashSetup;
namespace amrex { class Geometry; }

/**
 * @brief Fill an ADE plasma state from IMAS equilibrium/plasma data.
 *
 * This routine samples IMAS-derived plasma quantities on the AMReX grid and
 * fills the density, temperature, collision-frequency, and background magnetic
 * field arrays owned by `ADEPlasma`.
 *
 * @param setup IMAS sampling setup.
 * @param plasma Plasma model whose fields are initialized.
 * @param geom Geometry used to map grid cells to physical coordinates.
 */
void FillPlasmaFromImas(const ImashSetup& setup,
                        class ADEPlasma& plasma,
                        const amrex::Geometry& geom);
#endif

class EHJSolver;
class PMLSolver;

/**
 * @brief Auxiliary-differential-equation plasma response model coupled to an
 * electromagnetic FDTD solver.
 *
 * `ADEPlasma` evolves an electron momentum field and deposits the associated
 * plasma current into an attached `EHJSolver`. The update uses the local
 * electric field, prescribed/background magnetic field, electron density,
 * temperature, and collision frequency to compute the current response.
 *
 * The plasma model stores its fields on the same AMReX `BoxArray`,
 * `DistributionMapping`, and `Geometry` as the attached electromagnetic solver.
 * These grid objects are held by reference and must outlive the plasma object.
 *
 * The `sigma` field is also used as an activity/conductivity mask:
 *
 * \f[
 * \begin{array}{ll}
 * \sigma < 0 & \text{active plasma without Ohmic damping}, \\
 * \sigma = 0 & \text{inactive cell}, \\
 * \sigma > 0 & \text{active plasma with an added Ohmic current}.
 * \end{array}
 * \f]
 */
class ADEPlasma : public FieldOutputManager 
{
public:
    /**
     * @brief Construct an ADE plasma model attached to an electromagnetic
     * solver.
     *
     * The constructor allocates the electron momentum, density, temperature,
     * collision-frequency, conductivity-mask, and background magnetic-field
     * arrays. Plasma diagnostic fields are registered with the attached solver
     * so they are written with the electromagnetic output.
     *
     * @param solver Electromagnetic solver that supplies `E` and receives the
     * plasma current in `J`.
     */
    ADEPlasma(EHJSolver& solver);

    /**
     * @brief Return a read-only reference to the background magnetic field.
     *
     * @return Const reference to the three-component magnetic-field `MultiFab`.
     */
    amrex::MultiFab const& b() const { return B; }

    /**
     * @brief Return a mutable reference to the background magnetic field.
     *
     * @return Mutable reference to the three-component magnetic-field
     * `MultiFab`.
     */
    amrex::MultiFab& b() { return B; }

    /**
     * @brief Advance the electron momentum and accumulate plasma current.
     *
     * This method updates the electron momentum using the local electric field,
     * background magnetic field, collision frequency, density, and temperature.
     * The resulting current contribution is accumulated into the attached
     * solver's `J` field.
     *
     * The caller is responsible for ensuring that `solver.J` has the desired
     * initial value before this method is called. For a single plasma species,
     * this usually means zeroing `J`; for multiple current sources, accumulation
     * may be intentional.
     *
     * @param dt Timestep size.
     */
    void update(double dt);

    /**
     * @brief Fill the plasma conductivity mask from a PML scalar profile.
     *
     * The PML `sigmaJ` profile is copied into the ADE plasma `sigma` field and
     * scaled by the supplied maximum conductivity. Positive values activate an
     * additional Ohmic current contribution in the plasma update.
     *
     * @param pml PML solver providing the scalar `sigmaJ` profile.
     * @param sigma Maximum conductivity multiplier.
     */
    void fill_sigma(PMLSolver& pml, amrex::Real sigma = 50.0);

    /**
     * @brief Convert registered plasma diagnostic fields to dimensionless units.
     *
     * This scales density, temperature, collision frequency, and background
     * magnetic field by their normalization constants.
     */
    void remove_dimension();

    /**
     * @brief Restore physical units for registered plasma diagnostic fields.
     *
     * This reverses the scaling applied by `remove_dimension()`.
     */
    void add_dimension();

private:
    /**
     * @brief Electromagnetic solver coupled to this plasma response model.
     *
     * The plasma reads the solver electric field and accumulates current into
     * the solver current-density field.
     */
    EHJSolver& solver;

    /**
     * @brief AMReX box layout used by the plasma fields.
     *
     * Stored by reference; the referenced object must outlive the plasma model.
     */
    const amrex::BoxArray& BxA;

    /**
     * @brief Geometry descriptor associated with the plasma grid.
     *
     * Stored by reference; the referenced object must outlive the plasma model.
     */
    const amrex::Geometry& geom;

    /**
     * @brief Parallel distribution mapping for the plasma fields.
     *
     * Stored by reference; the referenced object must outlive the plasma model.
     */
    const amrex::DistributionMapping& dm;

    /**
     * @brief Electron momentum field.
     *
     * Components are ordered as `(px, py, pz)` using the same coordinate indices
     * as the electromagnetic fields.
     */
    amrex::MultiFab p;

protected:
    /**
     * @brief Electron number density.
     */
    amrex::MultiFab ne;

    /**
     * @brief Electron temperature.
     */
    amrex::MultiFab Te;

    /**
     * @brief Electron collision frequency.
     */
    amrex::MultiFab nu;

    /**
     * @brief Activity/conductivity mask used by the ADE update.
     *
     * The update interprets negative values as active plasma without Ohmic
     * damping, zero as inactive plasma, and positive values as active plasma
     * with an added Ohmic current.
     */
    amrex::MultiFab sigma;

    /**
     * @brief Prescribed/background magnetic flux density.
     *
     * Components are ordered as `(Bx, By, Bz)`.
     */
    amrex::MultiFab B;

    /**
     * @brief Maximum conductivity multiplier used when importing a PML profile.
     */
    amrex::Real sigma_max = 200.0;

#ifdef USE_IMAS
    friend void FillPlasmaFromImas(const ImashSetup& setup,
                                   ADEPlasma& plasma,
                                   const amrex::Geometry& geom);
#endif
};