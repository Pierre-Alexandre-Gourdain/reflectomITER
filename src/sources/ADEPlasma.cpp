#include "sources/ADEPlasma.H"
#include "utils/Stencils.H"
#include "fdtd/EHJSolver.H"
#include "fdtd/PMLSolver.H"

#include <cmath>

using namespace amrex;

ADEPlasma::ADEPlasma(EHJSolver& solver_)
    : solver(solver_),
      BxA(solver_.BxA),
      geom(solver_.geom),
      dm(solver_.dm)
{
    /*
     * The ADE model stores the electron momentum as a three-component field on
     * the same AMReX layout as the electromagnetic solver. The plasma
     * properties are scalar cell-centered fields, while B is a prescribed
     * three-component background magnetic field.
     */
    p.define(BxA, dm, 3, 0);

    ne.define(BxA, dm, 1, 0);
    Te.define(BxA, dm, 1, 0);
    nu.define(BxA, dm, 1, 0);
    sigma.define(BxA, dm, 1, 0);

    B.define(BxA, dm, 3, 0);

    /*
     * The plasma state is initialized to zero by default. sigma is initialized
     * to -1 so the update can distinguish ordinary plasma cells from optional
     * positive-conductivity PML or damping regions.
     */
    for (auto mf : {&p, &ne, &Te, &nu, &B})
    {
        mf->setVal(0.0);
    }

    sigma.setVal(-1.0);

    set_name("ADE_solver");

    /*
     * Register plasma diagnostics with the electromagnetic solver so field and
     * plasma quantities are written through the same output path.
     */
    solver.addField("ne", ne, geom, BxA, n_0, L_0);
    solver.addField("Te", Te, geom, BxA, T_0, L_0);
    solver.addField("nu", nu, geom, BxA, nu_0, L_0);
    solver.addField("B",  B,  geom, BxA, B_0, L_0);
}

void ADEPlasma::fill_sigma(PMLSolver& solver, Real sigma_)
{
    /*
     * Import the scalar PML profile into the plasma conductivity field. This
     * allows the ADE update to add an Ohmic current in regions where sigma is
     * positive, typically inside or near the absorbing layer.
     */
    sigma_max = sigma_;

    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto sigJ = solver.sigJ().const_array(mfi);
        auto sigma_arr = sigma.array(mfi);
        auto sig_max = sigma_max;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            sigma_arr(i,j,k) = sig_max * sigJ(i,j,k);
        });
    }
}

#if (false)
/*
 * Disabled alternative ADE update.
 *
 * This block preserves an older Strang-split/Boris-like momentum update:
 *
 *   1. half-step collisional damping,
 *   2. half-step electric acceleration,
 *   3. magnetic rotation,
 *   4. second half-step electric acceleration,
 *   5. final half-step collisional damping,
 *   6. current reconstruction from the updated momentum.
 *
 * It is kept as a reference implementation only. The active implementation
 * below uses an exact matrix-exponential ADE update for constant E, B, and nu
 * over one timestep.
 *
 * Do not enable this block without reviewing it first. In its current form it
 * contains legacy expressions and was not cleaned to the same standard as the
 * active update.
 */
void ADEPlasma::update(double dt)
{
    const double half_dt = 0.5 * dt;

    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto p_arr = p.array(mfi);
        auto E_arr = solver.E.const_array(mfi);
        auto J_arr = solver.J.array(mfi);

        auto ne_arr = ne.const_array(mfi);
        auto Te_arr = Te.const_array(mfi);
        auto B_arr = B.const_array(mfi);

        auto nu_arr = nu.const_array(mfi);
        auto sig_arr = sigma.const_array(mfi);

        AMREX_PARALLEL_FOR_3D(bx, i, j, k, {
            /*
             * The legacy convention uses sigma as an activity mask. Cells with
             * sigma exactly zero are skipped; negative sigma means active
             * plasma without Ohmic damping, and positive sigma means active
             * plasma with an additional Ohmic current.
             */
            if (!sig_arr(i,j,k)) return;

            const double nu_local = nu_arr(i,j,k);

            /*
             * This implementation assumes a collisional update. Collisionless
             * cells are skipped rather than handled by a separate limiting
             * expression.
             */
            if (nu_local < TINY) return;

            const double px_old = p_arr(i,j,k,X);
            const double py_old = p_arr(i,j,k,Y);
            const double pz_old = p_arr(i,j,k,Z);

            const double exp_half = exp(-nu_local * half_dt);

            /*
             * First half-step collisional damping.
             */
            const double p1x = exp_half * px_old;
            const double p1y = exp_half * py_old;
            const double p1z = exp_half * pz_old;

            /*
             * First half-step electric kick.
             */
            const double half_qdt = 0.5 * Q_e * dt;

            const double p2x = p1x + half_qdt * E_arr(i,j,k,X);
            const double p2y = p1y + half_qdt * E_arr(i,j,k,Y);
            const double p2z = p1z + half_qdt * E_arr(i,j,k,Z);

            /*
             * Magnetic rotation using the electron cyclotron vector.
             *
             * The legacy version uses the rest electron mass. The active update
             * below includes the temperature-dependent effective mass.
             */
            const double me = M_e;

            const double omX = Q_e * B_arr(i,j,k,X) / me;
            const double omY = Q_e * B_arr(i,j,k,Y) / me;
            const double omZ = Q_e * B_arr(i,j,k,Z) / me;

            const double om = sqrt(omX * omX + omY * omY + omZ * omZ);

            double p3x = p2x;
            double p3y = p2y;
            double p3z = p2z;

            if (om > TINY) {
                const double ux = omX / om;
                const double uy = omY / om;
                const double uz = omZ / om;

                const double theta = om * dt;
                const double cos_theta = cos(theta);
                const double sin_theta = sin(theta);

                /*
                 * Rodrigues rotation of the intermediate momentum around the
                 * magnetic-field direction. This is the central magnetic part
                 * of the split update.
                 */
                const double cross_x = uy * p2z - uz * p2y;
                const double cross_y = uz * p2x - ux * p2z;
                const double cross_z = ux * p2y - uy * p2x;

                const double udotp = ux * p2x + uy * p2y + uz * p2z;

                p3x = p2x * cos_theta
                    + cross_x * sin_theta
                    + ux * udotp * (1.0 - cos_theta);

                p3y = p2y * cos_theta
                    + cross_y * sin_theta
                    + uy * udotp * (1.0 - cos_theta);

                p3z = p2z * cos_theta
                    + cross_z * sin_theta
                    + uz * udotp * (1.0 - cos_theta);
            }

            /*
             * Second half-step electric kick.
             */
            const double p4x = p3x + half_qdt * E_arr(i,j,k,X);
            const double p4y = p3y + half_qdt * E_arr(i,j,k,Y);
            const double p4z = p3z + half_qdt * E_arr(i,j,k,Z);

            /*
             * Final half-step collisional damping and momentum storage.
             */
            p_arr(i,j,k,X) = exp_half * p4x;
            p_arr(i,j,k,Y) = exp_half * p4y;
            p_arr(i,j,k,Z) = exp_half * p4z;

            /*
             * Accumulate the plasma current associated with the updated
             * electron momentum.
             */
            J_arr(i,j,k,X) += Q_e * ne_arr(i,j,k) * p_arr(i,j,k,X) / me;
            J_arr(i,j,k,Y) += Q_e * ne_arr(i,j,k) * p_arr(i,j,k,Y) / me;
            J_arr(i,j,k,Z) += Q_e * ne_arr(i,j,k) * p_arr(i,j,k,Z) / me;

            /*
             * Positive sigma adds a local Ohmic current contribution.
             */
            if (sig_arr(i,j,k) > 0)
            {
                J_arr(i,j,k,X) += E_arr(i,j,k,X) * sig_arr(i,j,k);
                J_arr(i,j,k,Y) += E_arr(i,j,k,Y) * sig_arr(i,j,k);
                J_arr(i,j,k,Z) += E_arr(i,j,k,Z) * sig_arr(i,j,k);
            }
        });
    }
}
#else
void ADEPlasma::update(double dt)
{
    constexpr double small_nu = TINY;
    constexpr double small_om = TINY;

    for (MFIter mfi(BxA, dm, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();

        auto p_arr = p.array(mfi);
        auto E_arr = solver.E.const_array(mfi);
        auto J_arr = solver.J.array(mfi);

        auto ne_arr = ne.const_array(mfi);
        auto Te_arr = Te.const_array(mfi);
        auto B_arr  = B.const_array(mfi);
        auto nu_arr = nu.const_array(mfi);
        auto sig_arr = sigma.const_array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            /*
             * sigma acts as an activity mask as well as an optional Ohmic
             * conductivity. Cells with sigma exactly zero are skipped.
             */
            if (!sig_arr(i,j,k)) return;

            const double nu_local = nu_arr(i,j,k);
            if (nu_local < TINY) return;

            const double px0 = p_arr(i,j,k,X);
            const double py0 = p_arr(i,j,k,Y);
            const double pz0 = p_arr(i,j,k,Z);

            const double Ex = E_arr(i,j,k,X);
            const double Ey = E_arr(i,j,k,Y);
            const double Ez = E_arr(i,j,k,Z);

            /*
             * Effective electron mass used by the current model. The
             * temperature-dependent correction increases the inertia at high
             * electron temperature.
             */
            const double q_abs = std::abs(Q_e);
            const double m_el = M_e * std::sqrt(
                1.0 + 5.0 * q_abs * Te_arr(i,j,k) / (M_e * C_L2)
            );

            /*
             * Signed electron cyclotron vector. The charge sign is included in
             * Q_e, so the rotation direction follows the electron convention.
             */
            const double ox = Q_e * B_arr(i,j,k,X) / m_el;
            const double oy = Q_e * B_arr(i,j,k,Y) / m_el;
            const double oz = Q_e * B_arr(i,j,k,Z) / m_el;

            const double om2 = ox * ox + oy * oy + oz * oz;
            const double om  = std::sqrt(om2);

            double px_new = px0;
            double py_new = py0;
            double pz_new = pz0;

            if (om < small_om)
            {
                /*
                 * Unmagnetized ADE limit. The momentum equation reduces to a
                 * damped electric acceleration:
                 *
                 *     dp/dt = -nu p + q_e E.
                 */
                if (std::abs(nu_local) < small_nu)
                {
                    px_new = px0 + Q_e * dt * Ex;
                    py_new = py0 + Q_e * dt * Ey;
                    pz_new = pz0 + Q_e * dt * Ez;
                }
                else
                {
                    const double edamp = std::exp(-nu_local * dt);
                    const double c0 = (1.0 - edamp) / nu_local;

                    px_new = edamp * px0 + Q_e * c0 * Ex;
                    py_new = edamp * py0 + Q_e * c0 * Ey;
                    pz_new = edamp * pz0 + Q_e * c0 * Ez;
                }
            }
            else
            {
                /*
                 * Magnetized exact ADE update for constant E, B, and nu over
                 * the timestep. The operator K is defined by
                 *
                 *     K p = p x Omega,
                 *
                 * where Omega = q_e B / m_e.
                 */
                const double theta = om * dt;
                const double s = std::sin(theta);
                const double c = std::cos(theta);

                const double edamp = std::exp(-nu_local * dt);

                const double sin_over_om = s / om;
                const double one_minus_cos_over_om2 = (1.0 - c) / om2;

                const double Kpx = py0 * oz - pz0 * oy;
                const double Kpy = pz0 * ox - px0 * oz;
                const double Kpz = px0 * oy - py0 * ox;

                /*
                 * K^2 p is evaluated using the vector identity
                 *
                 *     K^2 p = Omega (Omega . p) - |Omega|^2 p.
                 */
                const double odotp = ox * px0 + oy * py0 + oz * pz0;

                const double K2px = ox * odotp - om2 * px0;
                const double K2py = oy * odotp - om2 * py0;
                const double K2pz = oz * odotp - om2 * pz0;

                /*
                 * Homogeneous momentum evolution:
                 *
                 *     exp[(-nu I + K) dt] p^n.
                 */
                double hx = px0
                          + sin_over_om * Kpx
                          + one_minus_cos_over_om2 * K2px;

                double hy = py0
                          + sin_over_om * Kpy
                          + one_minus_cos_over_om2 * K2py;

                double hz = pz0
                          + sin_over_om * Kpz
                          + one_minus_cos_over_om2 * K2pz;

                hx *= edamp;
                hy *= edamp;
                hz *= edamp;

                /*
                 * Exact electric forcing:
                 *
                 *     q_e integral_0^dt exp[(-nu I + K) tau] E d tau
                 *
                 * represented as q_e [c0 I + c1 K + c2 K^2] E.
                 */
                double c0;
                double c1;
                double c2;

                if (std::abs(nu_local) < small_nu)
                {
                    c0 = dt;
                    c1 = (1.0 - c) / om2;
                    c2 = (dt - s / om) / om2;
                }
                else
                {
                    const double den = nu_local * nu_local + om2;

                    c0 = (1.0 - edamp) / nu_local;

                    const double I_sin =
                        (om - edamp * (nu_local * s + om * c)) / den;

                    const double I_cos =
                        (nu_local + edamp * (-nu_local * c + om * s)) / den;

                    c1 = I_sin / om;
                    c2 = (c0 - I_cos) / om2;
                }

                const double KEx = Ey * oz - Ez * oy;
                const double KEy = Ez * ox - Ex * oz;
                const double KEz = Ex * oy - Ey * ox;

                const double odotE = ox * Ex + oy * Ey + oz * Ez;

                const double K2Ex = ox * odotE - om2 * Ex;
                const double K2Ey = oy * odotE - om2 * Ey;
                const double K2Ez = oz * odotE - om2 * Ez;

                px_new = hx + Q_e * (c0 * Ex + c1 * KEx + c2 * K2Ex);
                py_new = hy + Q_e * (c0 * Ey + c1 * KEy + c2 * K2Ey);
                pz_new = hz + Q_e * (c0 * Ez + c1 * KEz + c2 * K2Ez);
            }

            /*
             * Store the updated electron momentum and accumulate the
             * corresponding plasma current in the electromagnetic solver.
             */
            p_arr(i,j,k,X) = px_new;
            p_arr(i,j,k,Y) = py_new;
            p_arr(i,j,k,Z) = pz_new;

            J_arr(i,j,k,X) += Q_e * ne_arr(i,j,k) * px_new / m_el;
            J_arr(i,j,k,Y) += Q_e * ne_arr(i,j,k) * py_new / m_el;
            J_arr(i,j,k,Z) += Q_e * ne_arr(i,j,k) * pz_new / m_el;

            /*
             * A positive sigma adds a local Ohmic current contribution. This is
             * used to couple the ADE plasma response to damping or PML regions.
             */
            if (sig_arr(i,j,k) > 0)
            {
                J_arr(i,j,k,X) += E_arr(i,j,k,X) * sig_arr(i,j,k);
                J_arr(i,j,k,Y) += E_arr(i,j,k,Y) * sig_arr(i,j,k);
                J_arr(i,j,k,Z) += E_arr(i,j,k,Z) * sig_arr(i,j,k);
            }
        });
    }

    /*
     * The ADE update writes directly into solver.J. Fill its ghost cells before
     * any subsequent Maxwell update or divergence-cleaning operation reads J
     * across AMReX box boundaries.
     */
    solver.J.FillBoundary(geom.periodicity());
}
#endif

void ADEPlasma::remove_dimension()
{
    /*
     * Convert plasma diagnostic fields from physical units to the normalized
     * units used by output and post-processing.
     */
    FieldOutputManager::remove_dimension();

    ne.mult(1.0 / n_0, 0, ne.nComp(), ne.nGrow());
    Te.mult(1.0 / T_0, 0, Te.nComp(), Te.nGrow());
    nu.mult(1.0 / nu_0, 0, nu.nComp(), nu.nGrow());
    B .mult(1.0 / B_0, 0, B .nComp(), B .nGrow());
}

void ADEPlasma::add_dimension()
{
    /*
     * Restore plasma diagnostic fields to physical units after normalized
     * output processing.
     */
    FieldOutputManager::add_dimension();

    ne.mult(n_0, 0, ne.nComp(), ne.nGrow());
    Te.mult(T_0, 0, Te.nComp(), Te.nGrow());
    nu.mult(nu_0, 0, nu.nComp(), nu.nGrow());
    B .mult(B_0, 0, B .nComp(), B .nGrow());
}

#ifdef USE_IMAS
#include "io/ImasReader.H"

#include <AMReX_MFIter.H>
#include <AMReX_MultiFab.H>

namespace PlasmaCollisions
{
    inline double coulombLog_ei(double ne_m3, double Te_eV)
    {
        /*
         * Practical Coulomb logarithm estimate for hot, weakly coupled
         * plasmas. The density is converted to cm^-3 because the empirical
         * expression below is written in cgs density units.
         */
        if (ne_m3 <= 0.0 || Te_eV <= 0.0) {
            throw std::invalid_argument("ne_m3 and Te_eV must be > 0");
        }

        const double ne_cm3 = ne_m3 * 1.0e-6;

        double lnLambda = 23.0
                        - std::log(std::sqrt(ne_cm3) / std::pow(Te_eV, 1.5));

        /*
         * Clamp the estimate to a numerically reasonable range. This prevents
         * pathological IMAS samples from producing extreme collision rates.
         */
        lnLambda = std::clamp(lnLambda, 2.0, 30.0);

        return lnLambda;
    }

    inline double nu_ei(double ne_m3,
                        double Te_eV,
                        double Zeff,
                        double lnLambda)
    {
        /*
         * Spitzer-like practical formula for the electron-ion collision
         * frequency:
         *
         *     nu_ei [s^-1]
         *       = 2.91e-6 ne[cm^-3] Zeff lnLambda / Te[eV]^(3/2).
         */
        if (ne_m3 <= 0.0 || Te_eV <= 0.0 || Zeff <= 0.0 || lnLambda <= 0.0) {
            throw std::invalid_argument("All inputs must be > 0");
        }

        const double ne_cm3 = ne_m3 * 1.0e-6;
        return 2.91e-6 * ne_cm3 * Zeff * lnLambda / std::pow(Te_eV, 1.5);
    }

    inline double nu_ei(double ne_m3, double Te_eV, double Zeff = 1.0)
    {
        const double lnLambda = coulombLog_ei(ne_m3, Te_eV);
        return nu_ei(ne_m3, Te_eV, Zeff, lnLambda);
    }
}

using namespace PlasmaCollisions;

void FillPlasmaFromImas(const ImashSetup& setup,
                        ADEPlasma& plasma,
                        const amrex::Geometry& geom)
{
    using namespace amrex;

    /*
     * The IMAS initialization routine assumes scalar ne/Te/nu fields and a
     * three-component magnetic-field MultiFab.
     */
    AMREX_ALWAYS_ASSERT(plasma.ne.nComp() == 1);
    AMREX_ALWAYS_ASSERT(plasma.Te.nComp() == 1);
    AMREX_ALWAYS_ASSERT(plasma.nu.nComp() == 1);
    AMREX_ALWAYS_ASSERT(plasma.B.nComp()  == 3);

    /*
     * Sampling IMAS data is host-side work. Pinned host MultiFabs are used as
     * staging buffers before copying the initialized values back to the plasma
     * fields, which may live on device memory.
     */
    MFInfo info;
    info.SetArena(The_Pinned_Arena());

    MultiFab ne_h(plasma.ne.boxArray(), plasma.ne.DistributionMap(),
                  plasma.ne.nComp(), plasma.ne.nGrowVect(), info);
    MultiFab Te_h(plasma.Te.boxArray(), plasma.Te.DistributionMap(),
                  plasma.Te.nComp(), plasma.Te.nGrowVect(), info);
    MultiFab nu_h(plasma.nu.boxArray(), plasma.nu.DistributionMap(),
                  plasma.nu.nComp(), plasma.nu.nGrowVect(), info);
    MultiFab B_h (plasma.B .boxArray(), plasma.B .DistributionMap(),
                  plasma.B .nComp(), plasma.B .nGrowVect(), info);

    ne_h.setVal(0.0);
    Te_h.setVal(0.0);
    nu_h.setVal(0.0);
    B_h .setVal(0.0);

    const auto dx      = geom.CellSizeArray();
    const auto prob_lo = geom.ProbLoArray();

    for (MFIter mfi(ne_h, false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        auto ne_arr = ne_h.array(mfi);
        auto Te_arr = Te_h.array(mfi);
        auto nu_arr = nu_h.array(mfi);
        auto B_arr  = B_h.array(mfi);

        for (int k = bx.smallEnd(2); k <= bx.bigEnd(2); ++k) {
            for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {

                    /*
                     * Convert cell centers from normalized code coordinates to
                     * physical coordinates before sampling the IMAS data.
                     */
                    const double x = prob_lo[0] + (static_cast<double>(i) + 0.5) * dx[0];
                    const double y = prob_lo[1] + (static_cast<double>(j) + 0.5) * dx[1];
                    const double z = prob_lo[2] + (static_cast<double>(k) + 0.5) * dx[2];

                    const double x_phys = x * L_0;
                    const double y_phys = y * L_0;
                    const double z_phys = z * L_0;

                    const ImashSample s = SampleImashCartesian(setup,
                                                               x_phys,
                                                               y_phys,
                                                               z_phys);

                    ne_arr(i,j,k,0) = s.ne * DENSITY_FACTOR;
                    Te_arr(i,j,k,0) = s.Te;

                    if (s.ne * s.Te) {
                        nu_arr(i,j,k,0) = nu_ei(s.ne, s.Te);
                    }

                    B_arr(i,j,k,0) = s.Bx * B_FIELD_FACTOR;
                    B_arr(i,j,k,1) = s.By * B_FIELD_FACTOR;
                    B_arr(i,j,k,2) = s.Bz * B_FIELD_FACTOR;
                }
            }
        }
    }

    /*
     * Copy host-initialized IMAS data into the plasma MultiFabs used by the
     * solver.
     */
    plasma.ne.ParallelCopy(ne_h);
    plasma.Te.ParallelCopy(Te_h);
    plasma.nu.ParallelCopy(nu_h);
    plasma.B .ParallelCopy(B_h);
}
#endif