#include "fdtd/EHSolver.H"
#include "utils/Constants.H"
#include "utils/Stencils.H"

using namespace amrex;

EHSolver::EHSolver(const BoxArray& ba_in,
                   const Geometry& geom_in,
                   const DistributionMapping& dm_in,
                   const int nghost_)
    : BxA(ba_in),
      geom(geom_in),
      dm(dm_in),
      nghost(nghost_) 
{
    /*
     * The solver stores the electric and magnetic fields as cell-centered
     * three-component MultiFabs. The component order is assumed to match the
     * coordinate indices X, Y, and Z used throughout the stencil utilities.
     */
    E.define(BxA, dm, 3, nghost);
    H.define(BxA, dm, 3, nghost);

    set_name("FDTD_solver");

    /*
     * The fields start from a zero state. Physical initialization, source
     * deposition, or restart data loading should be applied by higher-level
     * code after construction.
     */
    E.setVal(0.0);
    H.setVal(0.0);

    /*
     * Register the fields with the common diagnostics/output infrastructure.
     * The normalization constants provide the physical scaling used when these
     * fields are written or inspected outside the solver.
     */
    addField("E", E, geom, BxA, E_0, L_0);
    addField("H", H, geom, BxA, H_0, L_0);
}

void EHSolver::updateE(Real dt)
{
    /*
     * The electric-field update uses the curl of H. Since the stencil reaches
     * into neighboring cells, the magnetic-field ghost cells must be valid
     * before entering the tile loop.
     */
    FillBoundaryH();

    dt *= eps0i;

    const auto ds = geom.CellSizeArray();
    const auto dxi = 1.0 / ds[X];
    const auto dyi = 1.0 / ds[Y];
    const auto dzi = 1.0 / ds[Z];

    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& ba = mfi.tilebox();

        auto ev = E.array(mfi);
        auto hv = H.const_array(mfi);

        AMREX_PARALLEL_FOR_3D(ba, i, j, k, {
            /*
             * Ampere's law in vacuum is applied in explicit form:
             *
             *     E^{n+1} = E^n + (dt / epsilon_0) curl(H).
             *
             * The backward-difference curl is paired with the forward
             * difference used in updateH(). This pairing is important for the
             * discrete FDTD operator and should not be changed independently.
             */
            ev(i,j,k,X) += dt * (bak_DY(hv,Z) - bak_DZ(hv,Y));
            ev(i,j,k,Y) += dt * (bak_DZ(hv,X) - bak_DX(hv,Z));
            ev(i,j,k,Z) += dt * (bak_DX(hv,Y) - bak_DY(hv,X));
        });
    }
}

void EHSolver::updateH(Real dt)
{
    /*
     * The magnetic-field update uses the curl of E. Electric-field ghost cells
     * must therefore be synchronized before evaluating the forward-difference
     * stencil at tile and box boundaries.
     */
    FillBoundaryE();

    dt *= mu0i;

    const auto ds = geom.CellSizeArray();
    const auto dxi = 1.0 / ds[X];
    const auto dyi = 1.0 / ds[Y];
    const auto dzi = 1.0 / ds[Z];

    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& ba = mfi.tilebox();

        auto hv = H.array(mfi);
        auto ev = E.const_array(mfi);

        AMREX_PARALLEL_FOR_3D(ba, i, j, k, {
            /*
             * Faraday's law is applied in explicit form:
             *
             *     H^{n+1} = H^n - (dt / mu_0) curl(E).
             *
             * This routine uses forward differences, complementing the
             * backward differences used in updateE().
             */
            hv(i,j,k,X) -= dt * (fwd_DY(ev,Z) - fwd_DZ(ev,Y));
            hv(i,j,k,Y) -= dt * (fwd_DZ(ev,X) - fwd_DX(ev,Z));
            hv(i,j,k,Z) -= dt * (fwd_DX(ev,Y) - fwd_DY(ev,X));
        });
    }
}

void EHSolver::FillBoundaryE()
{
    /*
     * AMReX exchanges periodic ghost cells and same-level neighboring data for
     * all electric-field components. Physical boundary conditions, if any, must
     * be handled elsewhere unless they are encoded through Geometry periodicity.
     */
    E.FillBoundary(geom.periodicity());
}

void EHSolver::FillBoundaryH()
{
    /*
     * AMReX exchanges periodic ghost cells and same-level neighboring data for
     * all magnetic-field components. This call is required before any curl
     * stencil that reads H across tile or box boundaries.
     */
    H.FillBoundary(geom.periodicity());
}