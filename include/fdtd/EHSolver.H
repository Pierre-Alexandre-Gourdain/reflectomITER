#ifndef EHSOLVER_H_
#define EHSOLVER_H_

#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include "io/FieldOutputManager.H"

class DivCleaning;

/**
 * @brief Vacuum finite-difference time-domain solver for the electric and
 * magnetic fields.
 *
 * `EHSolver` owns the cell-centered electromagnetic fields `E` and `H` and
 * advances them with an explicit FDTD update. The electric field is updated
 * from the discrete curl of `H`, and the magnetic field is updated from the
 * discrete curl of `E`.
 *
 * The fields are stored as three-component AMReX `MultiFab` objects using the
 * coordinate-component ordering `X`, `Y`, and `Z` defined in the solver
 * constants. The solver uses the provided AMReX `BoxArray`,
 * `DistributionMapping`, and `Geometry` by reference, so those objects must
 * outlive the solver.
 *
 * Although the vacuum solver does not use the current density directly, it
 * keeps a `J` member so derived and coupled solvers can share a common field
 * interface.
 */
class EHSolver : public FieldOutputManager {
public:
    /**
     * @brief Construct an electromagnetic FDTD solver on an AMReX grid.
     *
     * The constructor allocates the electric and magnetic field `MultiFab`s,
     * initializes them to zero, and registers them with the field-output
     * manager.
     *
     * @param ba AMReX box layout used by the solver.
     * @param geom Geometry descriptor associated with the grid.
     * @param dm Distribution mapping for parallel ownership of boxes.
     * @param nghost Number of ghost cells allocated for the field `MultiFab`s.
     */
    EHSolver(const amrex::BoxArray& ba,
             const amrex::Geometry& geom,
             const amrex::DistributionMapping& dm,
             const int nghost);

    /**
     * @brief Return a read-only reference to the current-density field.
     *
     * In the base vacuum solver this field is not used by the FDTD update, but
     * it is exposed for derived solvers and coupling components that require a
     * common field interface.
     *
     * @return Const reference to the current-density `MultiFab`.
     */
    amrex::MultiFab const& j() const { return J; }

    /**
     * @brief Return a mutable reference to the current-density field.
     *
     * @return Mutable reference to the current-density `MultiFab`.
     */
    amrex::MultiFab& j() { return J; }

    /**
     * @brief Return a read-only reference to the electric field.
     *
     * @return Const reference to the three-component electric-field `MultiFab`.
     */
    amrex::MultiFab const& e() const { return E; }

    /**
     * @brief Return a mutable reference to the electric field.
     *
     * @return Mutable reference to the three-component electric-field
     * `MultiFab`.
     */
    amrex::MultiFab& e() { return E; }

    /**
     * @brief Return a read-only reference to the magnetic field.
     *
     * @return Const reference to the three-component magnetic-field `MultiFab`.
     */
    amrex::MultiFab const& h() const { return H; }

    /**
     * @brief Return a mutable reference to the magnetic field.
     *
     * @return Mutable reference to the three-component magnetic-field
     * `MultiFab`.
     */
    amrex::MultiFab& h() { return H; }

    /**
     * @brief Destroy the solver.
     */
    virtual ~EHSolver() = default;

    /**
     * @brief Fill electric-field ghost cells.
     *
     * This synchronizes the electric-field `MultiFab` across AMReX box
     * boundaries using the geometry periodicity. Derived classes may override
     * this method to add physical boundary conditions or coupled-domain
     * exchanges.
     */
    virtual void FillBoundaryE();

    /**
     * @brief Fill magnetic-field ghost cells.
     *
     * This synchronizes the magnetic-field `MultiFab` across AMReX box
     * boundaries using the geometry periodicity. Derived classes may override
     * this method to add physical boundary conditions or coupled-domain
     * exchanges.
     */
    virtual void FillBoundaryH();

    /**
     * @brief Advance the electric field by one timestep.
     *
     * Applies the explicit FDTD electric-field update using the discrete curl
     * of the magnetic field:
     *
     * \f[
     *   \mathbf{E}^{n+1}
     *   =
     *   \mathbf{E}^{n}
     *   +
     *   \frac{\Delta t}{\epsilon_0}
     *   \nabla_h \times \mathbf{H}.
     * \f]
     *
     * @param dt Timestep size.
     */
    virtual void updateE(amrex::Real dt);

    /**
     * @brief Advance the magnetic field by one timestep.
     *
     * Applies the explicit FDTD magnetic-field update using the discrete curl
     * of the electric field:
     *
     * \f[
     *   \mathbf{H}^{n+1}
     *   =
     *   \mathbf{H}^{n}
     *   -
     *   \frac{\Delta t}{\mu_0}
     *   \nabla_h \times \mathbf{E}.
     * \f]
     *
     * @param dt Timestep size.
     */
    virtual void updateH(amrex::Real dt);

protected:
    /**
     * @brief AMReX box layout used by the solver.
     *
     * Stored by reference; the referenced object must outlive the solver.
     */
    const amrex::BoxArray& BxA;

    /**
     * @brief Geometry descriptor associated with the solver grid.
     *
     * Stored by reference; the referenced object must outlive the solver.
     */
    const amrex::Geometry& geom;

    /**
     * @brief Parallel distribution mapping for the solver boxes.
     *
     * Stored by reference; the referenced object must outlive the solver.
     */
    const amrex::DistributionMapping& dm;

    /**
     * @brief Three-component electric-field storage.
     *
     * Components are ordered as `(Ex, Ey, Ez)` using the coordinate indices
     * `X`, `Y`, and `Z`.
     */
    amrex::MultiFab E;

    /**
     * @brief Three-component magnetic-field storage.
     *
     * Components are ordered as `(Hx, Hy, Hz)` using the coordinate indices
     * `X`, `Y`, and `Z`.
     */
    amrex::MultiFab H;

    /**
     * @brief Three-component current-density storage.
     *
     * The base solver does not advance or use `J`, but the field is kept here
     * so derived solvers and coupling utilities can access current density
     * through a common interface.
     */
    amrex::MultiFab J;

    /**
     * @brief Number of ghost cells allocated for solver-owned fields.
     */
    const int nghost;

private:
    friend class ADEPlasma;
    friend class DivCleaning;
    friend class PointSource;
    friend class PMLSolver;
};

#endif