#ifndef DIVCLEANING_H_
#define DIVCLEANING_H_

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include "io/FieldOutputManager.H"

class EHJSolver;
class PMLSolver;

/**
 * @brief Hyperbolic divergence-cleaning operator for FDTD electromagnetic
 * fields.
 *
 * `DivCleaning` attaches to an existing `EHJSolver` and applies an auxiliary
 * scalar-potential correction to either the electric or magnetic field. When
 * cleaning the electric field, the target constraint is Gauss's law,
 *
 * \f[
 *   \nabla_h \cdot \mathbf{E} = \rho,
 * \f]
 *
 * where `rho` is advanced from the discrete continuity equation. When cleaning
 * the magnetic field, the target constraint is
 *
 * \f[
 *   \nabla_h \cdot \mathbf{H} = 0.
 * \f]
 *
 * The cleaner owns its auxiliary fields but operates directly on the `E`, `H`,
 * and `J` fields of the attached solver. The AMReX grid objects are stored by
 * reference and must outlive the cleaner.
 */
class DivCleaning : public FieldOutputManager {
public:
    /**
     * @brief Construct a divergence cleaner attached to an electromagnetic
     * solver.
     *
     * By default, the cleaner acts on the electric field and enforces Gauss's
     * law. Passing `clean_H = true` switches the cleaner to the magnetic
     * divergence constraint.
     *
     * @param solver_ Electromagnetic solver whose fields are cleaned.
     * @param clean_H If true, clean `div(H)` instead of the electric
     * Gauss-law residual.
     */
    DivCleaning(EHJSolver& solver_, bool clean_H = false);

    /**
     * @brief Fill ghost cells for the cleaner-owned auxiliary fields.
     *
     * For electric cleaning this synchronizes `rho`, `phi`, and `err`. For
     * magnetic cleaning this synchronizes `phi` and `err`.
     */
    void FillBoundary();

    /**
     * @brief Advance the charge density using the discrete continuity equation.
     *
     * This update is used only for electric-field cleaning. It computes
     * `rho_t = -div(J)` using the current-density field owned by the attached
     * solver.
     *
     * @param dt Timestep size.
     */
    void updateRho(amrex::Real dt);

    /**
     * @brief Advance the auxiliary cleaning potential.
     *
     * For electric cleaning, the potential is driven by the residual
     * `div(E) - rho`. For magnetic cleaning, it is driven by `div(H)`. The
     * potential is also damped on the cleaning timescale defined by the
     * implementation constants.
     *
     * @param dt Timestep size.
     */
    void updatePhi(amrex::Real dt);

    /**
     * @brief Apply the divergence-cleaning correction to the selected field.
     *
     * The correction subtracts the discrete gradient of the auxiliary potential
     * from either `E` or `H`, depending on whether this cleaner was constructed
     * for electric or magnetic divergence control.
     *
     * @param dt Timestep size.
     */
    void correct(amrex::Real dt);

    /**
     * @brief Apply the cleaning correction when this object is used as an
     * electric-field cleaner.
     *
     * This is a convenience wrapper around `correct()`.
     *
     * @param dt Timestep size.
     */
    void correctE(amrex::Real dt) { correct(dt); }

    /**
     * @brief Apply the cleaning correction when this object is used as a
     * magnetic-field cleaner.
     *
     * This is a convenience wrapper around `correct()`.
     *
     * @param dt Timestep size.
     */
    void correctH(amrex::Real dt) { correct(dt); }

    /**
     * @brief Fill electric-field ghost cells on the attached solver when needed.
     *
     * This method forwards to the attached solver only when the cleaner is
     * configured for electric-field cleaning.
     */
    void FillBoundaryE();

    /**
     * @brief Fill magnetic-field ghost cells on the attached solver when needed.
     *
     * This method forwards to the attached solver only when the cleaner is
     * configured for magnetic-field cleaning.
     */
    void FillBoundaryH();

private:
    /**
     * @brief Electromagnetic solver whose fields are cleaned.
     *
     * The cleaner does not operate on its own electromagnetic fields. It
     * modifies the fields stored in this attached solver.
     */
    EHJSolver& solver;

    /**
     * @brief AMReX box layout used by the cleaner.
     *
     * Stored by reference; the referenced object must outlive the cleaner.
     */
    const amrex::BoxArray& BxA;

    /**
     * @brief Geometry descriptor associated with the solver grid.
     *
     * Stored by reference; the referenced object must outlive the cleaner.
     */
    const amrex::Geometry& geom;

    /**
     * @brief Parallel distribution mapping for the cleaner fields.
     *
     * Stored by reference; the referenced object must outlive the cleaner.
     */
    const amrex::DistributionMapping& dm;

    /**
     * @brief Number of ghost cells allocated for cleaner-owned fields.
     */
    const int nghost;

    /**
     * @brief Charge-density estimate used by electric-field cleaning.
     *
     * This field is allocated only when the cleaner enforces Gauss's law for
     * the electric field.
     */
    amrex::MultiFab rho;

    /**
     * @brief Auxiliary scalar potential used to propagate and damp divergence
     * error.
     */
    amrex::MultiFab phi;

    /**
     * @brief Logarithmic diagnostic of the divergence-cleaning residual.
     *
     * The implementation stores `log10(abs(residual) + TINY)` for visualization
     * and diagnostics.
     */
    amrex::MultiFab err;

    /**
     * @brief Selects electric or magnetic divergence cleaning.
     *
     * If true, the cleaner targets `div(E) - rho`. If false, it targets
     * `div(H)`.
     */
    bool clean_E;
};

#endif