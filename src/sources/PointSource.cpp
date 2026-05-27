#include "sources/PointSource.H"

#include <cmath>

#include <AMReX_MFIter.H>
#include <AMReX_MultiFabUtil.H>

using namespace amrex;

namespace
{
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real source_profile_value(Real r,
                              Real width,
                              Source::ProfileType profile,
                              Real supergaussian_order) noexcept
    {
        /*
         * A nonpositive width degenerates the source to a single-cell impulse.
         * The radial profile is then nonzero only exactly at the source center.
         */
        if (width <= Real(0.0)) {
            return (r == Real(0.0)) ? Real(1.0) : Real(0.0);
        }

        switch (profile)
        {
            case Source::ProfileType::Gaussian:
            {
                /*
                 * Smooth Gaussian profile with width interpreted in physical
                 * units after the anisotropic coordinate scaling has been
                 * applied.
                 */
                Real x = r / width;
                return std::exp(-Real(0.5) * x * x);
            }

            case Source::ProfileType::SuperGaussian:
            {
                /*
                 * Super-Gaussian profile. Orders smaller than one are clamped
                 * to one so the profile remains well behaved.
                 */
                Real m = (supergaussian_order < Real(1.0)) ? Real(1.0)
                                                           : supergaussian_order;
                Real x = r / width;
                return std::exp(-Real(0.5) * std::pow(x, Real(2.0) * m));
            }

            case Source::ProfileType::Heaviside:
            {
                /*
                 * Compact top-hat profile. The width is used as the physical
                 * radius of the support.
                 */
                return (r <= width) ? Real(1.0) : Real(0.0);
            }
        }

        return Real(0.0);
    }
}

PointSource::PointSource(EHSolver& solver_,
                         SourceType type_,
                         Real x_, Real y_, Real z_,
                         Real amplitude_,
                         Real frequency_,
                         Real phase_,
                         Real width_,
                         ProfileType profile_,
                         Real supergaussian_order_,
                         Vector<Real> const& direction_,
                         const Geometry& geom_,
                         Vector<Real> const& shape_,
                         amrex::Real t_on,
                         amrex::Real t_off,
                         amrex::Real t_rise)

    : EHsolver(&solver_),
      type(type_),
      x0(x_), y0(y_), z0(z_),
      amplitude(amplitude_),
      frequency(frequency_),
      phase(phase_),
      width(width_),
      profile(profile_),
      supergaussian_order(supergaussian_order_),
      geom(geom_)
{
    /*
     * Cache the physical lower corner and grid spacing. Source deposition later
     * converts between physical coordinates and cell indices without repeatedly
     * querying the AMReX Geometry object.
     */
    auto const plo = geom_.ProbLoArray();
    auto const ddx = geom_.CellSizeArray();

    for (int d = 0; d < 3; ++d) {
        prob_lo[d] = plo[d];
        dx[d]      = ddx[d];
    }

    /*
     * Normalize the requested source polarization direction. A zero direction
     * vector is treated as a zero vector after normalization protection, so it
     * deposits no component-wise field unless the caller supplies a direction.
     */
    Real norm = 0.0;
    for (int d = 0; d < 3; ++d) {
        norm += direction_[d] * direction_[d];
    }

    norm = std::sqrt(norm);

    if (norm <= 0.0) norm = 1.0;

    for (int d = 0; d < 3; ++d) {
        dir[d] = direction_[d] / norm;
    }

    /*
     * The shape factors stretch the radial coordinate independently in each
     * direction. TINY prevents division by zero when a caller requests a
     * vanishing extent in one direction.
     */
    for (int d = 0; d < 3; ++d) {
        shape[d] = std::abs(shape_[d]) + TINY;
    }

    /*
     * Store the source turn-on and turn-off interval in increasing time order.
     * A value of -1 for either endpoint disables the temporal envelope.
     */
    t_start = std::min(t_on, t_off);
    t_stop  = std::max(t_on, t_off);
    t_switch = t_rise;

    /*
     * For oscillatory sources, prevent the switching time from becoming much
     * shorter than the wave period. This avoids an abrupt envelope that would
     * inject unnecessary high-frequency content into the simulation.
     */
	if (frequency > Real(0.0)) {
		Real min_switch = Real(1.0) / frequency / Real(5.0);
		if (t_switch < min_switch) {
			t_switch = min_switch;
		}
	}
}

void PointSource::update(Real time, Real dt)
{
    Real val = amplitude;

    if (frequency > Real(0.0)) {
        Real omega = Real(2.0) * Real(M_PI) * frequency;
        Real arg = omega * time + phase * Real(M_PI / 180.0);

        /*
         * A current source is deposited as an oscillating current density. A
         * direct E or H source is deposited as an increment to the field, so
         * the sinusoid is differentiated in time and multiplied by dt.
         */
        if (type == SourceType::J) {
            val *= std::sin(arg);
        } else {
            val *= omega * dt * std::cos(arg);
        }

        /*
         * Optional smooth source envelope. The product of two tanh ramps turns
         * the source on near t_start and off near t_stop without a discontinuity
         * in the deposited signal.
         */
        if (t_start != Real(-1.0) && t_stop != Real(-1.0)) {
            Real envelope_on =
                (std::tanh((time - t_start) / t_switch) + Real(1.0)) / Real(2.0);

            Real envelope_off =
                (std::tanh((t_stop - time) / t_switch) + Real(1.0)) / Real(2.0);

            val *= envelope_on * envelope_off;
        }
    }

    /*
     * Deposit into the selected field. Each field is checked before use so the
     * same source object can be attached safely to solvers that may not allocate
     * all possible source targets.
     */
    switch (type)
    {
        case SourceType::J:
            if (EHsolver && EHsolver->J.ok()) setPoint(EHsolver->J, val);
            break;

        case SourceType::E:
            if (EHsolver && EHsolver->E.ok()) setPoint(EHsolver->E, val);
            break;

        case SourceType::H:
            if (EHsolver && EHsolver->H.ok()) setPoint(EHsolver->H, val);
            break;

        default:
            break;
    }
}

void PointSource::setPoint(MultiFab& mf, Real val)
{
    /*
     * The profile width is interpreted in physical units. Smooth profiles are
     * truncated at three widths to keep the deposition local while retaining
     * the dominant part of the Gaussian or super-Gaussian support.
     */
    Real rmax = (profile == ProfileType::Heaviside)
                    ? width
                    : Real(3.0) * width;

    if (rmax < Real(0.0)) rmax = Real(0.0);

    /*
     * Convert the physical support radius to an index-space search box. The
     * anisotropic shape factors can stretch the support differently along each
     * coordinate direction.
     */
    int nx = static_cast<int>(std::ceil(rmax / dx[0] * std::max(Real(1.0), shape[0])));
    int ny = static_cast<int>(std::ceil(rmax / dx[1] * std::max(Real(1.0), shape[1])));
    int nz = static_cast<int>(std::ceil(rmax / dx[2] * std::max(Real(1.0), shape[2])));

    /*
     * Locate the cell containing the source center. The source box is then
     * clipped against each AMReX tile before deposition.
     */
    int ic = static_cast<int>(std::floor((x0 - prob_lo[0]) / dx[0]));
    int jc = static_cast<int>(std::floor((y0 - prob_lo[1]) / dx[1]));
    int kc = static_cast<int>(std::floor((z0 - prob_lo[2]) / dx[2]));

    IntVect lo(ic - nx, jc - ny, kc - nz);
    IntVect hi(ic + nx, jc + ny, kc + nz);
    Box srcBox(lo, hi);

    for (MFIter mfi(mf, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        auto arr = mf.array(mfi);

        auto const x00 = x0;
        auto const y00 = y0;
        auto const z00 = z0;
        auto const plo = prob_lo;
        auto const ddx = dx;
        auto const w   = width;
        auto const p   = profile;
        auto const sg  = supergaussian_order;
        auto const d   = dir;
        auto const sh  = shape;

        Box overlap = bx & srcBox;

        if (overlap.ok())
        {
            AMREX_PARALLEL_FOR_3D(overlap, i, j, k,
            {
                /*
                 * Evaluate the source profile at the cell center. The radial
                 * coordinate is computed after applying the anisotropic shape
                 * scaling, so the support can represent elongated or flattened
                 * source regions.
                 */
                Real x = plo[0] + (Real(i) + Real(0.5)) * ddx[0];
                Real y = plo[1] + (Real(j) + Real(0.5)) * ddx[1];
                Real z = plo[2] + (Real(k) + Real(0.5)) * ddx[2];

                Real rx = (x - x00) / sh[0];
                Real ry = (y - y00) / sh[1];
                Real rz = (z - z00) / sh[2];
                Real r  = std::sqrt(rx * rx + ry * ry + rz * rz);

                Real profile_value = source_profile_value(r, w, p, sg);

                /*
                 * Deposit the source into all vector components according to
                 * the normalized polarization direction.
                 */
                arr(i,j,k,0) += val * profile_value * d[0];
                arr(i,j,k,1) += val * profile_value * d[1];
                arr(i,j,k,2) += val * profile_value * d[2];
            });
        }
    }
}