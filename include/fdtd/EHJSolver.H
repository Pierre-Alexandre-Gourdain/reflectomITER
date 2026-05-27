#ifndef EHJSOLVER_H_
#define EHJSOLVER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Array4.H>

#include "fdtd/EHSolver.H"

/**
 * @brief Electromagnetic FDTD solver with an explicit current-density source.
 *
 * `EHJSolver` extends `EHSolver` by allocating and advancing with a
 * three-component current-density field `J`. The magnetic-field update is
 * inherited from `EHSolver`, while the electric-field update adds the current
 * source term in Ampere's law:
 *
 * \f[
 *   \mathbf{E}^{n+1}
 *   =
 *   \mathbf{E}^{n}
 *   +
 *   \frac{\Delta t}{\epsilon_0}
 *   \left(
 *     \nabla_h \times \mathbf{H}
 *     -
 *     \mathbf{J}
 *   \right).
 * \f]
 *
 * The field arrays use the same AMReX grid layout, distribution mapping, and
 * geometry as the base electromagnetic solver. The current density is stored
 * with the same component ordering as the electric and magnetic fields.
 */
class EHJSolver : public EHSolver {
public:
    /**
     * @brief Construct an electromagnetic solver with current-density storage.
     *
     * The constructor initializes the base `EHSolver`, allocates the
     * three-component current-density `MultiFab`, initializes it to zero, and
     * registers it with the field-output manager.
     *
     * @param ba AMReX box layout used by the solver.
     * @param geom Geometry descriptor associated with the grid.
     * @param dm Distribution mapping for parallel ownership of boxes.
     * @param nghost Number of ghost cells allocated for the field `MultiFab`s.
     */
    EHJSolver(const amrex::BoxArray& ba,
              const amrex::Geometry& geom,
              const amrex::DistributionMapping& dm,
              const int nghost);

    /**
     * @brief Advance the electric field by one timestep including current.
     *
     * This override first applies the vacuum electric-field update from
     * `EHSolver`, then applies the explicit current-density source term
     * \f$-\Delta t\,\mathbf{J}/\epsilon_0\f$.
     *
     * The current-density ghost cells are filled before the update so that the
     * source term is valid across AMReX box and tile boundaries.
     *
     * @param dt Timestep size.
     */
    void updateE(amrex::Real dt) override;

    /**
     * @brief Fill current-density ghost cells.
     *
     * This synchronizes the current-density `MultiFab` across AMReX box
     * boundaries using the geometry periodicity. Derived classes may override
     * the electric or magnetic boundary-fill logic separately if they implement
     * coupled-domain exchanges or physical boundary conditions.
     */
    void FillBoundaryJ();
};

#endif