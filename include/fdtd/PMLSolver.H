#ifndef PML_SOLVER_H_
#define PML_SOLVER_H_

#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Array4.H>

#include "io/FieldOutputManager.H"
#include "utils/HelperFunctions.H"
#include "fdtd/EHJSolver.H"

#include <cmath>

/**
 * @brief Convolutional perfectly matched layer solver coupled to an EHJ FDTD
 * solver.
 *
 * `PMLSolver` implements an absorbing-layer update for the electromagnetic
 * fields using CPML-style auxiliary memory variables. The PML owns its own
 * `E`, `H`, and `J` fields through `EHJSolver`, but it is coupled to a
 * separate physical-domain `EHJSolver` through explicit two-way boundary
 * exchanges.
 *
 * The solver stores directional damping profiles for electric and magnetic
 * updates, auxiliary `psi` memory fields for split derivative terms, and
 * correction fields `phi_E` and `phi_H` that are applied after the standard
 * FDTD update.
 *
 * A zero-thickness PML is treated as inactive. This allows driver code to keep
 * the same update sequence regardless of whether the absorbing layer is
 * enabled.
 */
class PMLSolver : public EHJSolver {
public:
    /**
     * @brief Construct a PML solver coupled to a physical electromagnetic
     * solver.
     *
     * The PML uses its own AMReX grid layout, geometry, and distribution
     * mapping. The attached solver supplies the physical fields with which the
     * PML exchanges boundary data.
     *
     * @param solver Physical-domain EHJ solver coupled to this PML.
     * @param ba AMReX box layout for the PML region.
     * @param geom Geometry descriptor for the PML grid.
     * @param dm Distribution mapping for parallel ownership of PML boxes.
     * @param thickness Number of cells used for the absorbing layer. A value of
     * zero disables the PML.
     * @param sigma Maximum damping strength used to build the PML profile.
     */
    PMLSolver(EHJSolver& solver,
              const amrex::BoxArray& ba,
              const amrex::Geometry& geom,
              const amrex::DistributionMapping& dm,
              int thickness = 10,
              amrex::Real sigma = 100.0);

    /**
     * @brief Return a read-only reference to the electric damping profile.
     *
     * @return Const reference to the three-component `sigmaE` field.
     */
    amrex::MultiFab const& sigE() const { return sigmaE; }

    /**
     * @brief Return a mutable reference to the electric damping profile.
     *
     * @return Mutable reference to the three-component `sigmaE` field.
     */
    amrex::MultiFab& sigE() { return sigmaE; }

    /**
     * @brief Return a read-only reference to the magnetic damping profile.
     *
     * @return Const reference to the three-component `sigmaH` field.
     */
    amrex::MultiFab const& sigH() const { return sigmaH; }

    /**
     * @brief Return a mutable reference to the magnetic damping profile.
     *
     * @return Mutable reference to the three-component `sigmaH` field.
     */
    amrex::MultiFab& sigH() { return sigmaH; }

    /**
     * @brief Return a read-only reference to the scalar current/PML profile.
     *
     * @return Const reference to the one-component `sigmaJ` field.
     */
    amrex::MultiFab const& sigJ() const { return sigmaJ; }

    /**
     * @brief Return a mutable reference to the scalar current/PML profile.
     *
     * @return Mutable reference to the one-component `sigmaJ` field.
     */
    amrex::MultiFab& sigJ() { return sigmaJ; }

    /**
     * @brief Exchange electric-field data between the physical solver and the
     * PML solver.
     *
     * The physical solver and the PML solver own separate `E` fields on
     * different AMReX layouts. This routine performs the required two-way
     * `ParallelCopy` so both domains see consistent data in their coupling
     * region.
     */
    void FullBoundaryExchangeE();

    /**
     * @brief Exchange magnetic-field data between the physical solver and the
     * PML solver.
     *
     * The physical solver and the PML solver own separate `H` fields on
     * different AMReX layouts. This routine performs the required two-way
     * `ParallelCopy` so both domains see consistent data in their coupling
     * region.
     */
    void FullBoundaryExchangeH();

    /**
     * @brief Exchange current-density data between the physical solver and the
     * PML solver.
     *
     * This keeps the current source term consistent across the physical-domain
     * and absorbing-layer coupling region.
     */
    void FullBoundaryExchangeJ();

    /**
     * @brief Update the electric CPML auxiliary correction field.
     *
     * This routine updates the electric-side CPML memory variables using
     * derivatives of `H`, then assembles `phi_E`. The correction field is later
     * applied by `updateE()`.
     *
     * @param dt Timestep size.
     */
    void updatePhiE(amrex::Real dt);

    /**
     * @brief Update the magnetic CPML auxiliary correction field.
     *
     * This routine updates the magnetic-side CPML memory variables using
     * derivatives of `E`, then assembles `phi_H`. The correction field is later
     * applied by `updateH()`.
     *
     * @param dt Timestep size.
     */
    void updatePhiH(amrex::Real dt);

    /**
     * @brief Advance the PML electric field by one timestep.
     *
     * The update first synchronizes magnetic-field data with the coupled
     * physical solver, then applies the standard EHJ electric-field update,
     * followed by the CPML electric correction.
     *
     * @param dt Timestep size.
     */
    void updateE(amrex::Real dt) override;

    /**
     * @brief Advance the PML magnetic field by one timestep.
     *
     * The update first synchronizes electric-field data with the coupled
     * physical solver, then applies the standard magnetic-field update,
     * followed by the CPML magnetic correction.
     *
     * @param dt Timestep size.
     */
    void updateH(amrex::Real dt) override;

    /**
     * @brief Allocate and initialize the PML damping and auxiliary fields.
     *
     * This builds the directional damping profiles, initializes the CPML memory
     * coefficients, and zeroes the auxiliary fields used by the absorbing-layer
     * update.
     */
    void initializeSigma();

    /**
     * @brief Convert registered PML diagnostic fields to dimensionless units.
     *
     * This overrides the output scaling behavior for the PML damping profiles,
     * which have inverse-length dimensions.
     */
    void remove_dimension();

    /**
     * @brief Restore physical units for registered PML diagnostic fields.
     *
     * This reverses the scaling applied by `remove_dimension()`.
     */
    void add_dimension();

private:
    /**
     * @brief Physical-domain solver coupled to this PML.
     *
     * The PML exchanges field data with this solver but stores and updates its
     * own PML fields separately.
     */
    EHJSolver& solver;

    /**
     * @brief Number of cells in the absorbing layer.
     *
     * A value of zero disables the PML update.
     */
    int thickness;

    /**
     * @brief Maximum damping coefficient used in the PML profile.
     */
    amrex::Real sigma_max;

    /**
     * @brief Polynomial order used to shape the PML damping profile.
     */
    int poly_order = 2;

    /**
     * @brief Electric CPML correction field applied after the standard
     * electric-field update.
     */
    amrex::MultiFab phi_E;

    /**
     * @brief Magnetic CPML correction field applied after the standard
     * magnetic-field update.
     */
    amrex::MultiFab phi_H;

    /**
     * @brief CPML memory variables associated with the electric-field update
     * and the Ex component.
     */
    amrex::MultiFab psi_Ex;

    /**
     * @brief CPML memory variables associated with the electric-field update
     * and the Ey component.
     */
    amrex::MultiFab psi_Ey;

    /**
     * @brief CPML memory variables associated with the electric-field update
     * and the Ez component.
     */
    amrex::MultiFab psi_Ez;

    /**
     * @brief CPML memory variables associated with the magnetic-field update
     * and the Hx component.
     */
    amrex::MultiFab psi_Hx;

    /**
     * @brief CPML memory variables associated with the magnetic-field update
     * and the Hy component.
     */
    amrex::MultiFab psi_Hy;

    /**
     * @brief CPML memory variables associated with the magnetic-field update
     * and the Hz component.
     */
    amrex::MultiFab psi_Hz;

    /**
     * @brief Electric CPML exponential decay coefficient.
     */
    amrex::MultiFab bE;

    /**
     * @brief Electric CPML derivative accumulation coefficient.
     */
    amrex::MultiFab cE;

    /**
     * @brief Magnetic CPML exponential decay coefficient.
     */
    amrex::MultiFab bH;

    /**
     * @brief Magnetic CPML derivative accumulation coefficient.
     */
    amrex::MultiFab cH;

    /**
     * @brief Directional electric damping profile.
     *
     * Components correspond to damping in the `x`, `y`, and `z` directions.
     */
    amrex::MultiFab sigmaE;

    /**
     * @brief Directional magnetic damping profile.
     *
     * Components correspond to damping in the `x`, `y`, and `z` directions.
     */
    amrex::MultiFab sigmaH;

    /**
     * @brief Scalar PML profile associated with current-density handling or
     * diagnostics.
     */
    amrex::MultiFab sigmaJ;

    /**
     * @brief Small CPML stability parameter added to the damping coefficient.
     *
     * This corresponds to the CPML alpha parameter. The default value disables
     * the additional alpha damping.
     */
    amrex::Real a_pml = 0.0;
};

/**
 * @brief Compute the exponential CPML decay coefficient.
 *
 * This returns
 *
 * \f[
 *   b = \exp(-\gamma \Delta t).
 * \f]
 *
 * @param gamma Effective CPML damping coefficient.
 * @param dt Timestep size.
 * @return Exponential decay coefficient.
 */
AMREX_GPU_HOST_DEVICE
inline amrex::Real safe_b_from_gamma_dt(amrex::Real gamma,
                                        amrex::Real dt) noexcept
{
    return std::exp(-gamma * dt);
}

/**
 * @brief Compute the CPML derivative accumulation coefficient.
 *
 * This returns
 *
 * \f[
 *   c = \frac{1 - \exp(-\gamma \Delta t)}{\gamma}.
 * \f]
 *
 * For nonpositive `gamma`, the limiting value `dt` is returned.
 *
 * @param gamma Effective CPML damping coefficient.
 * @param dt Timestep size.
 * @return Derivative accumulation coefficient.
 */
AMREX_GPU_HOST_DEVICE
inline amrex::Real safe_c_from_gamma_dt(amrex::Real gamma,
                                        amrex::Real dt) noexcept
{
    if (gamma <= 0.0) return dt;

#if defined(__CUDA_ARCH__)
    amrex::Real e = std::exp(-gamma * dt);
    return (1.0 - e) / gamma;
#else
    amrex::Real x = -gamma * dt;
    amrex::Real num = -std::expm1(x);
    return num / gamma;
#endif
}

#endif