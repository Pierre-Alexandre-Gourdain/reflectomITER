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
#include "utils/Stencils.H"
#include "utils/HelperFunctions.H"
#include "fdtd/PMLSolver.H"
#include "fdtd/EHJSolver.H"

#include <cmath>

using namespace amrex;

PMLSolver::PMLSolver(EHJSolver& solver_,
                     const BoxArray& ba_,
                     const Geometry& geom_,
                     const DistributionMapping& dm_,
                     int thickness_,
                     Real sigma_)
    : EHJSolver(ba_, geom_, dm_, solver_.nghost),
      solver(solver_),
      thickness(thickness_),
      sigma_max(sigma_)
{
    /*
     * A zero-thickness PML is treated as an inactive object. This lets the
     * higher-level driver keep the same call structure whether or not an
     * absorbing layer is present.
     */
    if (!thickness) return;

    initializeSigma();

    set_name("pml");

    /*
     * Match the output directory of the physical solver so that PML
     * diagnostics are written together with the main electromagnetic fields.
     */
    output_dir = solver.output_dir;

    /*
     * The sigma fields are useful diagnostics: they show where the absorbing
     * layer is active and how strongly each coordinate direction is damped.
     */
    addField("sigmaE", sigmaE, geom, BxA);
    addField("sigmaH", sigmaH, geom, BxA);
    addField("sigmaJ", sigmaJ, geom, BxA);
}

void PMLSolver::remove_dimension()
{
    /*
     * FieldOutputManager removes the dimensions from the registered fields.
     * The PML damping profiles are inverse-length quantities, so they require
     * their own scaling when switching to dimensionless output.
     */
    FieldOutputManager::remove_dimension();

    sigmaE.mult(L_0, 0, sigmaE.nComp(), sigmaE.nGrow());
    sigmaH.mult(L_0, 0, sigmaH.nComp(), sigmaH.nGrow());
    sigmaJ.mult(L_0, 0, sigmaJ.nComp(), sigmaJ.nGrow());
}

void PMLSolver::add_dimension()
{
    /*
     * Restore the physical inverse-length scaling of the PML damping profiles
     * after dimensionless output or diagnostic processing.
     */
    FieldOutputManager::add_dimension();

    sigmaE.mult(1.0 / L_0, 0, sigmaE.nComp(), sigmaE.nGrow());
    sigmaH.mult(1.0 / L_0, 0, sigmaH.nComp(), sigmaH.nGrow());
    sigmaJ.mult(1.0 / L_0, 0, sigmaJ.nComp(), sigmaJ.nGrow());
}

void PMLSolver::initializeSigma()
{
    if (!thickness) return;

    /*
     * The CPML auxiliary fields store the history terms associated with split
     * derivatives. Each vector component has its own memory contribution so the
     * PML can damp waves according to the coordinate direction in which they
     * enter the absorbing layer.
     */
    phi_E.define(BxA, dm, 3, nghost);
    phi_H.define(BxA, dm, 3, nghost);
    phi_E.setVal(0.0);
    phi_H.setVal(0.0);

    psi_Ex.define(BxA, dm, 3, nghost);
    psi_Ey.define(BxA, dm, 3, nghost);
    psi_Ez.define(BxA, dm, 3, nghost);

    psi_Hx.define(BxA, dm, 3, nghost);
    psi_Hy.define(BxA, dm, 3, nghost);
    psi_Hz.define(BxA, dm, 3, nghost);

    psi_Ex.setVal(0.0);
    psi_Ey.setVal(0.0);
    psi_Ez.setVal(0.0);

    psi_Hx.setVal(0.0);
    psi_Hy.setVal(0.0);
    psi_Hz.setVal(0.0);

    /*
     * bE/cE and bH/cH are the time-integration coefficients for the CPML memory
     * variables. They are recomputed from sigma and dt before each update.
     */
    bE.define(BxA, dm, 3, nghost);
    cE.define(BxA, dm, 3, nghost);
    bH.define(BxA, dm, 3, nghost);
    cH.define(BxA, dm, 3, nghost);

    bE.setVal(1.0);
    cE.setVal(0.0);
    bH.setVal(1.0);
    cH.setVal(0.0);

    /*
     * sigmaE and sigmaH store directional damping profiles for the electric and
     * magnetic CPML updates. sigmaJ is scalar because it is intended to damp or
     * identify current-density behavior in the PML region.
     */
    sigmaE.define(BxA, dm, 3, nghost);
    sigmaH.define(BxA, dm, 3, nghost);
    sigmaJ.define(BxA, dm, 1, nghost);

    sigmaE.setVal(0.0);
    sigmaH.setVal(0.0);
    sigmaJ.setVal(0.0);

    /*
     * The damping profile is built relative to the physical solver domain. The
     * PML BoxArray can extend outside this domain; cells farther into the layer
     * receive larger damping according to the selected polynomial order.
     */
    const Box domain = solver.BxA.minimalBox();
    const Real sig_max = sigma_max;
    const int th = thickness;
    const auto order = poly_order;

    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto sE = sigmaE.array(mfi);
        auto sH = sigmaH.array(mfi);
        auto sJ = sigmaJ.array(mfi);

        /*
         * The two passes assign the directional E/H damping first and then
         * construct a scalar PML indicator for J in cells not already covered
         * by the directional sigma profile.
         */
        for (int m = 1; m >= 0; m--)
        {
            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                Real dx = 0.0;
                Real dy = 0.0;
                Real dz = 0.0;

                if (i < domain.smallEnd()[X] - th / 2 * m) {
                    dx = domain.smallEnd()[X] - th / 2 * m - i;
                } else if (i > domain.bigEnd()[X] + th / 2 * m) {
                    dx = i - domain.bigEnd()[X] - th / 2 * m;
                }

                if (j < domain.smallEnd()[Y] - th / 2 * m) {
                    dy = domain.smallEnd()[Y] - th / 2 * m - j;
                } else if (j > domain.bigEnd()[Y] + th / 2 * m) {
                    dy = j - domain.bigEnd()[Y] - th / 2 * m;
                }

                if (k < domain.smallEnd()[Z] - th / 2 * m) {
                    dz = domain.smallEnd()[Z] - th / 2 * m - k;
                } else if (k > domain.bigEnd()[Z] + th / 2 * m) {
                    dz = k - domain.bigEnd()[Z] - th / 2 * m;
                }

                Real sx = (dx > 0.0) ? pow(dx / th * 2, order) : 0.0;
                Real sy = (dy > 0.0) ? pow(dy / th * 2, order) : 0.0;
                Real sz = (dz > 0.0) ? pow(dz / th * 2, order) : 0.0;

                if (m == 1)
                {
                    sE(i,j,k,X) = sig_max * sx;
                    sE(i,j,k,Y) = sig_max * sy;
                    sE(i,j,k,Z) = sig_max * sz;

                    sH(i,j,k,X) = sE(i,j,k,X);
                    sH(i,j,k,Y) = sE(i,j,k,Y);
                    sH(i,j,k,Z) = sE(i,j,k,Z);
                }
                else
                {
                    /*
                     * sigmaJ is only populated where the directional sigma
                     * vector is effectively zero. This makes it a scalar marker
                     * for the inner PML transition region.
                     */
                    if (sE(i,j,k,X) * sE(i,j,k,X)
                      + sE(i,j,k,Y) * sE(i,j,k,Y)
                      + sE(i,j,k,Z) * sE(i,j,k,Z) < TINY)
                    {
                        sJ(i,j,k,X) = std::sqrt(sx * sx + sy * sy + sz * sz)
                                    + TINY;
                    }
                }
            });
        }
    }
}

void PMLSolver::updatePhiE(Real dt)
{
    if (!thickness) return;

    /*
     * phi_E is driven by derivatives of H, so H ghost cells must be valid
     * before evaluating the backward-difference curl contributions.
     */
    FillBoundaryH();

    /*
     * Recompute the CPML memory coefficients for the current timestep. This
     * allows dt to vary between calls without invalidating the PML update.
     */
    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto sE = sigmaE.const_array(mfi);
        auto bEarr = bE.array(mfi);
        auto cEarr = cE.array(mfi);
        auto a = a_pml;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            for (int d = 0; d < 3; ++d) {
                Real sigma_val = sE(i,j,k,d);
                Real gamma = sigma_val + a;

                Real b = safe_b_from_gamma_dt(gamma, dt);
                Real c = safe_c_from_gamma_dt(gamma, dt);

                bEarr(i,j,k,d) = b;
                cEarr(i,j,k,d) = sigma_val * c;
            }
        });
    }

    const auto ds = geom.CellSizeArray();
    const auto dxi = 1.0 / ds[X];
    const auto dyi = 1.0 / ds[Y];
    const auto dzi = 1.0 / ds[Z];
	
    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto Hv = H.const_array(mfi);

        auto psiEx = psi_Ex.array(mfi);
        auto psiEy = psi_Ey.array(mfi);
        auto psiEz = psi_Ez.array(mfi);
        auto phiEarr = phi_E.array(mfi);

        auto bEarr = bE.const_array(mfi);
        auto cEarr = cE.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            /*
             * The CPML memory variables decay exponentially and then accumulate
             * the derivative terms that would otherwise appear directly in the
             * Maxwell curl operator.
             */
            for (int d = 0; d < 3; ++d) {
                Real bval = bEarr(i,j,k,d);

                psiEx(i,j,k,d) *= bval;
                psiEy(i,j,k,d) *= bval;
                psiEz(i,j,k,d) *= bval;
            }

            psiEx(i,j,k,Y) += bak_DY(Hv, Z) * cEarr(i,j,k,Y);
            psiEx(i,j,k,Z) += bak_DZ(Hv, Y) * cEarr(i,j,k,Z);

            psiEy(i,j,k,Z) += bak_DZ(Hv, X) * cEarr(i,j,k,Z);
            psiEy(i,j,k,X) += bak_DX(Hv, Z) * cEarr(i,j,k,X);

            psiEz(i,j,k,X) += bak_DX(Hv, Y) * cEarr(i,j,k,X);
            psiEz(i,j,k,Y) += bak_DY(Hv, X) * cEarr(i,j,k,Y);

            /*
             * Assemble the PML correction using the same sign pattern as the
             * electric-field curl update. updateE() applies this correction
             * after the usual EHJ Maxwell step.
             */
            phiEarr(i,j,k,X) = psiEx(i,j,k,Y) - psiEx(i,j,k,Z);
            phiEarr(i,j,k,Y) = psiEy(i,j,k,Z) - psiEy(i,j,k,X);
            phiEarr(i,j,k,Z) = psiEz(i,j,k,X) - psiEz(i,j,k,Y);
        });
    }
}

void PMLSolver::updateE(Real dt)
{
    if (!thickness) return;

    /*
     * Synchronize H between the physical solver and the PML solver before the
     * electric update. The PML update uses the local PML field storage but must
     * remain consistent with the adjacent physical domain.
     */
    FullBoundaryExchangeH();

    /*
     * Apply the standard EHJ electric-field advance first:
     *
     *     E <- E + dt/epsilon_0 curl(H) - dt/epsilon_0 J.
     *
     * The PML correction is then applied as an additional source term.
     */
    EHJSolver::updateE(dt);

    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto Earr = E.array(mfi);
        auto phiEarr = phi_E.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Earr(i,j,k,X) -= dt * phiEarr(i,j,k,X);
            Earr(i,j,k,Y) -= dt * phiEarr(i,j,k,Y);
            Earr(i,j,k,Z) -= dt * phiEarr(i,j,k,Z);
        });
    }

    /*
     * Push the updated PML electric field back to the coupled solver so the
     * neighboring physical domain uses the corrected boundary data.
     */
    FullBoundaryExchangeE();
}

void PMLSolver::updatePhiH(Real dt)
{
    if (!thickness) return;

    /*
     * phi_H is driven by derivatives of E, so E ghost cells must be valid
     * before evaluating the forward-difference curl contributions.
     */
    FillBoundaryE();

    /*
     * Recompute the magnetic CPML memory coefficients for the current timestep.
     */
    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto sH = sigmaH.const_array(mfi);
        auto bHarr = bH.array(mfi);
        auto cHarr = cH.array(mfi);
        auto a = a_pml;

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            for (int d = 0; d < 3; ++d) {
                Real sigma_val = sH(i,j,k,d);
                Real gamma = sigma_val + a;

                Real b = safe_b_from_gamma_dt(gamma, dt);
                Real c = safe_c_from_gamma_dt(gamma, dt);

                bHarr(i,j,k,d) = b;
                cHarr(i,j,k,d) = sigma_val * c;
            }
        });
    }

    const auto ds = geom.CellSizeArray();
    const auto dxi = 1.0 / ds[X];
    const auto dyi = 1.0 / ds[Y];
    const auto dzi = 1.0 / ds[Z];

    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto Ev = E.const_array(mfi);

        auto psiHx = psi_Hx.array(mfi);
        auto psiHy = psi_Hy.array(mfi);
        auto psiHz = psi_Hz.array(mfi);
        auto phiHarr = phi_H.array(mfi);

        auto bHarr = bH.const_array(mfi);
        auto cHarr = cH.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            /*
             * The magnetic CPML memory variables mirror the electric memory
             * update, but use the forward-difference stencil associated with
             * the H update.
             */
            for (int d = 0; d < 3; ++d) {
                Real bval = bHarr(i,j,k,d);

                psiHx(i,j,k,d) *= bval;
                psiHy(i,j,k,d) *= bval;
                psiHz(i,j,k,d) *= bval;
            }

            psiHx(i,j,k,Y) += fwd_DY(Ev, Z) * cHarr(i,j,k,Y);
            psiHx(i,j,k,Z) += fwd_DZ(Ev, Y) * cHarr(i,j,k,Z);

            psiHy(i,j,k,Z) += fwd_DZ(Ev, X) * cHarr(i,j,k,Z);
            psiHy(i,j,k,X) += fwd_DX(Ev, Z) * cHarr(i,j,k,X);

            psiHz(i,j,k,X) += fwd_DX(Ev, Y) * cHarr(i,j,k,X);
            psiHz(i,j,k,Y) += fwd_DY(Ev, X) * cHarr(i,j,k,Y);

            /*
             * Assemble the PML correction using the sign pattern appropriate
             * for the magnetic-field curl update.
             */
            phiHarr(i,j,k,X) = psiHx(i,j,k,Y) - psiHx(i,j,k,Z);
            phiHarr(i,j,k,Y) = psiHy(i,j,k,Z) - psiHy(i,j,k,X);
            phiHarr(i,j,k,Z) = psiHz(i,j,k,X) - psiHz(i,j,k,Y);
        });
    }
}

void PMLSolver::updateH(Real dt)
{
    if (!thickness) return;

    /*
     * Synchronize E between the physical solver and the PML solver before the
     * magnetic update. H depends on the curl of E, so stale E boundary data
     * would contaminate the absorbing-layer update.
     */
    FullBoundaryExchangeE();

    /*
     * Apply the standard magnetic-field advance before adding the PML
     * correction.
     */
    EHSolver::updateH(dt);

    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();

        auto Hv = H.array(mfi);
        auto phiHarr = phi_H.const_array(mfi);

        ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            Hv(i,j,k,X) += dt * phiHarr(i,j,k,X);
            Hv(i,j,k,Y) += dt * phiHarr(i,j,k,Y);
            Hv(i,j,k,Z) += dt * phiHarr(i,j,k,Z);
        });
    }

    /*
     * Push the updated PML magnetic field back to the coupled solver so the
     * neighboring physical domain receives the corrected boundary data.
     */
    FullBoundaryExchangeH();
}

void PMLSolver::FullBoundaryExchangeH()
{
    if (!thickness) return;

    /*
     * The PML owns its own H MultiFab but is coupled to the H field in the
     * physical solver. This exchange keeps overlapping regions and ghost cells
     * synchronized between the two objects.
     */
    FillBoundaryH();

    for (auto [src, dst] : {
        std::pair(&solver.H, &H),
        std::pair(&H, &solver.H),
    }) {
        dst->ParallelCopy(*src,
                          0,
                          0,
                          dst->nComp(),
                          IntVect(0),
                          dst->nGrowVect(),
                          geom.periodicity());
    }

    FillBoundaryH();
}

void PMLSolver::FullBoundaryExchangeE()
{
    if (!thickness) return;

    /*
     * The PML owns its own E MultiFab but is coupled to the E field in the
     * physical solver. This exchange keeps overlapping regions and ghost cells
     * synchronized between the two objects.
     */
    FillBoundaryE();

    for (auto [src, dst] : {
        std::pair(&solver.E, &E),
        std::pair(&E, &solver.E),
    }) {
        dst->ParallelCopy(*src,
                          0,
                          0,
                          dst->nComp(),
                          IntVect(0),
                          dst->nGrowVect(),
                          geom.periodicity());
    }

    FillBoundaryE();
}

void PMLSolver::FullBoundaryExchangeJ()
{
    if (!thickness) return;

    /*
     * Current density is exchanged in the same way as E and H so the PML and
     * physical solver see a consistent source term near the absorbing layer.
     */
    FillBoundaryJ();

    for (auto [src, dst] : {
        std::pair(&solver.J, &J),
        std::pair(&J, &solver.J),
    }) {
        dst->ParallelCopy(*src,
                          0,
                          0,
                          dst->nComp(),
                          IntVect(0),
                          dst->nGrowVect(),
                          geom.periodicity());
    }

    FillBoundaryJ();
}