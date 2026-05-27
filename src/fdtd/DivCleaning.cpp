#include "utils/Constants.H"
#include "utils/Stencils.H"
#include "fdtd/EHJSolver.H"
#include "fdtd/DivCleaning.H"

#include <AMReX_MultiFab.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Geometry.H>

using namespace amrex;

DivCleaning::DivCleaning(EHJSolver& solver_, bool clean_H)
    :BxA(solver_.BxA),
     geom(solver_.geom),
     dm(solver_.dm),
     nghost(solver_.nghost), 
	 solver(solver_)
{
    /*
     * The same cleaner class is used for either electric or magnetic
     * divergence control. Passing clean_H = true selects div(H) cleaning;
     * otherwise the cleaner enforces Gauss's law for E.
     */
	
    clean_E = !clean_H;

    if (clean_E)
    {
        /*
         * Electric-field cleaning requires a charge-density estimate because
         * the target constraint is div(E) = rho. The auxiliary potential phi
         * propagates and damps the constraint error, while err stores a
         * diagnostic measure of the residual.
         */
        for (auto mf : {&rho, &phi, &err})
        {
            mf->define(BxA, dm, 1, nghost);
            mf->setVal(0.0);
        }

        /*
         * The charge density and logarithmic electric-divergence error are
         * registered on the parent solver so they appear with the physical
         * fields in diagnostics and plotfile output.
         */
        solver.addField("rho", rho, geom, BxA, rho_0, L_0);
        solver.addField("Eerr", err, geom, BxA, 1, 1);
    }
    else
    {
        /*
         * Magnetic cleaning does not require a charge-density field. The
         * magnetic constraint is div(H) = 0, so only the auxiliary potential
         * and residual diagnostic are needed.
         */
        for (auto mf : {&phi, &err})
        {
            mf->define(BxA, dm, 1, nghost);
            mf->setVal(0.0);
        }

        /*
         * The divergence error is stored on a log10 scale and is therefore
         * treated as dimensionless in the output registry.
         */
        solver.addField("Herr", err, geom, BxA, 1, 1);
    }

    set_name("div_cleaner");
}

void DivCleaning::updateRho(Real dt)
{
    /*
     * Charge density is evolved from the discrete continuity equation:
     *
     *     rho_t = - div(J).
     *
     * If the attached solver has no valid current-density field, there is no
     * source from which to update rho.
     */
    if (!solver.J.ok()) return;

    /*
     * Important: the current density used below is solver.J, not the inherited
     * J field from EHJSolver. Therefore the boundary fill must be applied to
     * the attached solver's current field.
     */
    solver.FillBoundaryJ();

    const auto ds = geom.CellSizeArray();
    const auto dxi = 1.0 / ds[X];
    const auto dyi = 1.0 / ds[Y];
    const auto dzi = 1.0 / ds[Z];
	
    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto J_arr = solver.J.const_array(mfi);
        auto rho_arr = rho.array(mfi);

        AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
            /*
             * The same backward-divergence convention used for div(E) is used
             * here so that the charge update is discretely consistent with the
             * electric Gauss-law residual.
             */
            Real divJ = bak_DX(J_arr, X)
                      + bak_DY(J_arr, Y)
                      + bak_DZ(J_arr, Z);

            rho_arr(i, j, k) -= dt * divJ;
        });
    }
}

void DivCleaning::updatePhi(Real dt)
{
    const auto ds = geom.CellSizeArray();
    const auto dxi = 1.0 / ds[X];
    const auto dyi = 1.0 / ds[Y];
    const auto dzi = 1.0 / ds[Z];

    if (clean_E)
    {
        /*
         * Before evaluating the electric constraint, update rho from the
         * current density and synchronize both the physical E field and the
         * cleaner's scalar fields.
         */
        updateRho(dt);
        FillBoundaryE();
        FillBoundary();

        for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            auto E_arr = solver.E.const_array(mfi);
            auto phi_arr = phi.array(mfi);
            auto rho_arr = rho.const_array(mfi);
            auto err_arr = err.array(mfi);

            AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
                /*
                 * The electric residual is the discrete Gauss-law error. The
                 * auxiliary potential is driven by this residual and damped on
                 * the timescale Tau.
                 */
                Real divE = bak_DX(E_arr, X)
                          + bak_DY(E_arr, Y)
                          + bak_DZ(E_arr, Z);

                Real err_val = divE - rho_arr(i, j, k);

                /*
                 * Store the residual on a logarithmic scale for visualization.
                 * TINY prevents log10(0) while preserving small-error trends.
                 */
                err_arr(i, j, k) = std::log10(std::abs(err_val) + TINY);

                phi_arr(i, j, k) -= dt * (err_val * C_h * C_h
                                        + phi_arr(i, j, k) / Tau);
            });
        }
    }
    else
    {
        /*
         * Magnetic cleaning targets div(H) = 0. The H ghost cells and cleaner
         * scalar fields must be valid before the residual is evaluated.
         */
        FillBoundaryH();
        FillBoundary();

        for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            auto H_arr = solver.H.const_array(mfi);
            auto phi_arr = phi.array(mfi);
            auto err_arr = err.array(mfi);

            AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
                /*
                 * The magnetic residual uses the forward-divergence convention.
                 * This should remain paired with the gradient convention used
                 * in correct() so that the cleaner acts as a consistent
                 * hyperbolic/parabolic divergence-control operator.
                 */
                Real divH = fwd_DX(H_arr, X)
                          + fwd_DY(H_arr, Y)
                          + fwd_DZ(H_arr, Z);

                Real err_val = divH;

                err_arr(i, j, k) = std::log10(std::abs(err_val) + TINY);

                phi_arr(i, j, k) -= dt * (err_val * C_h * C_h
                                        + phi_arr(i, j, k) / Tau);
            });
        }
    }
}

void DivCleaning::correct(Real dt)
{
    /*
     * The correction step applies the gradient of the auxiliary potential to
     * the field being cleaned. Phi ghost cells must be synchronized because the
     * gradient stencil reads neighboring cells.
     */
    FillBoundary();

    const auto ds = geom.CellSizeArray();
    const auto dxi = 1.0 / ds[X];
    const auto dyi = 1.0 / ds[Y];
    const auto dzi = 1.0 / ds[Z];

    if (clean_E)
    {
        for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            auto E_arr = solver.E.array(mfi);
            auto phi_arr = phi.const_array(mfi);

            AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
                /*
                 * The electric correction removes the gradient part associated
                 * with the propagated Gauss-law residual:
                 *
                 *     E <- E - dt grad(phi).
                 */
                E_arr(i, j, k, X) -= dt * fwd_dx(phi_arr);
                E_arr(i, j, k, Y) -= dt * fwd_dy(phi_arr);
                E_arr(i, j, k, Z) -= dt * fwd_dz(phi_arr);
            });
        }
    }
    else
    {
        for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();

            auto H_arr = solver.H.array(mfi);
            auto phi_arr = phi.const_array(mfi);

            AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
                /*
                 * The magnetic correction applies the auxiliary-potential
                 * gradient to H. Each field component must use the derivative
                 * in the matching coordinate direction.
                 */
                H_arr(i, j, k, X) -= dt * bak_dx(phi_arr);
                H_arr(i, j, k, Y) -= dt * bak_dy(phi_arr);
                H_arr(i, j, k, Z) -= dt * bak_dz(phi_arr);
            });
        }
    }
}

void DivCleaning::FillBoundary()
{
    if (clean_E)
    {
        /*
         * Electric cleaning owns rho, phi, and err. All three are synchronized
         * together so diagnostics and subsequent stencil operations see
         * consistent ghost-cell data.
         */
        for (auto* mf : {&rho, &phi, &err})
        {
            mf->FillBoundary(geom.periodicity());
        }
    }
    else
    {
        /*
         * Magnetic cleaning owns only phi and err. Physical open-boundary
         * treatment, if required, should be added here in addition to same-level
         * AMReX ghost exchange.
         */
        for (auto* mf : {&phi, &err})
        {
            mf->FillBoundary(geom.periodicity());
        }
    }
}

void DivCleaning::FillBoundaryE()
{
    /*
     * DivCleaning acts on the fields owned by the attached solver rather than
     * on the E/H storage inherited from EHJSolver. Forward boundary-fill calls
     * to the attached solver to keep the operated-on field valid.
     */
    if (clean_E)
    {
        solver.FillBoundaryE();
    }
}

void DivCleaning::FillBoundaryH()
{
    /*
     * DivCleaning acts on the fields owned by the attached solver rather than
     * on the E/H storage inherited from EHJSolver. Forward boundary-fill calls
     * to the attached solver to keep the operated-on field valid.
     */
    if (!clean_E)
    {
        solver.FillBoundaryH();
    }
}