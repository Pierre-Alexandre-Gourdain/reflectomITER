#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>
#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_Array4.H>

#include "utils/Constants.H"
#include "utils/HelperFunctions.H"
#include "fdtd/EHJSolver.H"

using namespace amrex;

EHJSolver::EHJSolver(const BoxArray& ba_in,
                     const Geometry& geom_in,
                     const DistributionMapping& dm_in,
                     const int nghost)
    : EHSolver(ba_in, geom_in, dm_in, nghost)
{
    /*
     * EHJSolver extends the vacuum E/H solver with a three-component current
     * density field. The current is stored on the same BoxArray and
     * DistributionMapping as E and H so that the field update can be performed
     * locally inside each AMReX tile.
     */
    J.define(BxA, dm, 3, nghost);
    J.setVal(0.0);

    /*
     * Register J with the common field/output infrastructure. The scaling
     * constants are used by diagnostics and plotfile output to convert the
     * stored code units into the expected physical normalization.
     */
    addField("J", J, geom, BxA, J_0, L_0);
}

void EHJSolver::updateE(Real dt)
{
    /*
     * The current-density ghost cells are filled before the electric-field
     * update so that J is valid at tile and box boundaries. This matters when
     * J is produced by another solver component or communicated across AMReX
     * boxes.
     */
    FillBoundaryJ();

    /*
     * First apply the vacuum Ampere update implemented by the base class:
     *
     *     E <- E + (dt / epsilon_0) curl(H).
     *
     * The current contribution is applied below as a separate explicit source
     * term.
     */
    EHSolver::updateE(dt);

    dt *= eps0i;

    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& ba = mfi.tilebox();

        auto ev = E.array(mfi);
        auto jv = J.const_array(mfi);

        AMREX_PARALLEL_FOR_3D(ba, i, j, k, {
            /*
             * Apply the current-density term in Ampere's law:
             *
             *     E <- E - (dt / epsilon_0) J.
             *
             * The update is component-wise because J is already stored in the
             * same component ordering as E.
             */
            ev(i,j,k,X) -= dt * jv(i,j,k,X);
            ev(i,j,k,Y) -= dt * jv(i,j,k,Y);
            ev(i,j,k,Z) -= dt * jv(i,j,k,Z);
        });
    }
}

void EHJSolver::FillBoundaryJ()
{
    /*
     * Synchronize current-density ghost cells using the geometry periodicity.
     * Non-periodic physical boundary conditions, if required, should be applied
     * by the code that owns the current model or boundary-condition logic.
     */
    J.FillBoundary(geom.periodicity());
}