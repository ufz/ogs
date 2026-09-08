// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "DriftFluxModel.h"

#include <spdlog/fmt/fmt.h>

#include <algorithm>
#include <cmath>
#include <string>

#include "MaterialLib/PhysicalConstant.h"
#include "NumLib/Exceptions.h"

namespace MaterialPropertyLib
{
/// Slip parameter \f$S = C_0 (x + (1 - x) \rho_v / \rho_l)\f$ of the
/// drift-flux closure, see computeVapourVoidFraction().
static double closureSlipParameter(DriftFluxState const& state)
{
    return state.C_0 *
           (state.dryness + (1 - state.dryness) * state.vapour_water_density /
                                state.liquid_water_density);
}

/// The same state with the drift aligned with the mixture flow, which is the
/// closure computeVapourVoidFraction() actually solves.
static DriftFluxState alignDriftWithFlow(DriftFluxState const& state)
{
    DriftFluxState aligned = state;
    aligned.u_gu = alignedDriftFluxVelocity(state.u_gu, state.v_mix);
    return aligned;
}

double alignedDriftFluxVelocity(double const u_gu, double const v_mix)
{
    return std::copysign(u_gu, v_mix);
}

double driftFluxProfileParameter(double const dryness)
{
    return 1 + 0.12 * (1 - dryness);
}

VoidFractionQuadratic voidFractionQuadratic(DriftFluxState const& state)
{
    double const S = closureSlipParameter(state);
    double const delta =
        state.vapour_water_density - state.liquid_water_density;

    return {
        -S * delta * state.v_mix,
        state.v_mix * (state.dryness * delta - S * state.liquid_water_density) -
            state.vapour_water_density * state.u_gu,
        state.dryness * state.v_mix * state.liquid_water_density};
}

double voidFractionResidual(double const alpha, DriftFluxState const& state)
{
    double const rho_mix = alpha * state.vapour_water_density +
                           (1 - alpha) * state.liquid_water_density;
    double const S = closureSlipParameter(state);

    return state.liquid_water_density *
           (rho_mix * state.v_mix * (state.dryness - alpha * S) -
            alpha * state.vapour_water_density * state.u_gu);
}

std::string voidFractionClosureDiagnostics(DriftFluxState const& state)
{
    DriftFluxState const aligned = alignDriftWithFlow(state);

    // computeVapourVoidFraction() accepts a root only from [0, alpha_max], so
    // the sign change has to be reported over that interval. Over [0, 1] it
    // says nothing: the second root of the closure is always larger than one,
    // and the residual changes sign between alpha_max and one whenever the
    // admissible interval holds no root at all.
    double const alpha_max =
        std::min(1., aligned.dryness / closureSlipParameter(aligned));

    return fmt::format(
        "dryness {:g}, liquid density {:g} kg/m^3, vapour density {:g} kg/m^3, "
        "profile parameter {:g}, aligned drift flux velocity {:g} m/s. The "
        "closure residual is {:g} at a void fraction of zero and {:g} at the "
        "upper bound {:g} of the admissible interval; a root exists in between "
        "only if the two have opposite signs.",
        aligned.dryness, aligned.liquid_water_density,
        aligned.vapour_water_density, aligned.C_0, aligned.u_gu,
        voidFractionResidual(0., aligned),
        voidFractionResidual(alpha_max, aligned), alpha_max);
}

double driftFluxVelocity(double const dryness, double const temperature,
                         double const vapour_water_density,
                         double const liquid_water_density)
{
    // The rounded value the original implementation of this closure used, not
    // standard gravity, whose defined value is 9.80665 m/s^2. It stays rounded
    // because sharpening it would move every existing two-phase result, and it
    // stays local because a rounded stand-in is not the physical constant and
    // so does not belong in PhysicalConstant.h beside the critical point.
    constexpr double gravity = 9.81;  // m/s^2

    // Both caps are written as comparisons rather than with std::max, which
    // returns its first argument for a NaN second one and would hand out a
    // finite drift flux velocity for a non-finite state. The NaN is kept so
    // that computeVapourVoidFraction() sees it and aborts the assembly.
    double const temperature_ratio =
        1 - temperature /
                MaterialLib::PhysicalConstant::CriticalPoint::TemperatureWater;
    double const reduced_temperature =
        temperature_ratio < 0 ? 0. : temperature_ratio;
    double const sigma_gl = 0.2358 * std::pow(reduced_temperature, 1.256) *
                            (1 - 0.625 * reduced_temperature);

    double const buoyancy =
        gravity * sigma_gl * (liquid_water_density - vapour_water_density);
    double const drift_buoyancy = buoyancy < 0 ? 0. : buoyancy;

    return 1.18 * (1 - dryness) * std::pow(drift_buoyancy, 0.25) /
           std::pow(liquid_water_density, 0.5);
}

DriftFluxState driftFluxState(double const dryness, double const temperature,
                              double const vapour_water_density,
                              double const liquid_water_density,
                              double const v_mix)
{
    return {.dryness = dryness,
            .vapour_water_density = vapour_water_density,
            .liquid_water_density = liquid_water_density,
            .v_mix = v_mix,
            .C_0 = driftFluxProfileParameter(dryness),
            .u_gu = alignedDriftFluxVelocity(
                driftFluxVelocity(dryness, temperature, vapour_water_density,
                                  liquid_water_density),
                v_mix)};
}

std::optional<double> computeVapourVoidFraction(DriftFluxState const& state)
{
    auto const& [dryness, vapour_water_density, liquid_water_density, v_mix,
                 C_0, u_gu] = state;

    // Single phase states are exact, no closure is needed. They are settled
    // before the state validation below because the density of the absent
    // phase does not enter the result: it is the one evaluated far off the
    // saturation line and hence the one that may be non-positive, and a pure
    // liquid or pure vapour section must not abort the assembly for it. A NaN
    // dryness satisfies neither comparison and reaches the validation.
    if (dryness <= 0)
    {
        return 0.;
    }
    if (dryness >= 1)
    {
        return 1.;
    }

    // All arguments are computed from the current solution iterate, so a
    // non-physical value is a diverging global Newton step, not a broken
    // input. Aborting the assembly ends the nonlinear iteration and lets the
    // time stepping repeat the step under either solver, whereas OGS_FATAL
    // would end the run; see the documentation of this function.
    //
    // The phase densities are evaluated on the IAPWS-IF97 region 4 saturation
    // line, which is extrapolated outside of its pressure range 611.213 Pa to
    // 22.064 MPa and then may return non-positive or NaN values. The negated
    // comparison is deliberate; `density <= 0` lets a NaN density pass.
    if (!(liquid_water_density > 0) || !(vapour_water_density > 0))
    {
        throw NumLib::AssemblyException(fmt::format(
            "Non-positive phase density in the vapour void fraction closure: "
            "liquid density {:g} kg/m^3, vapour density {:g} kg/m^3.",
            liquid_water_density, vapour_water_density));
    }
    if (!std::isfinite(dryness) || !std::isfinite(v_mix) ||
        !std::isfinite(u_gu) || !std::isfinite(C_0))
    {
        throw NumLib::AssemblyException(fmt::format(
            "Non-finite state in the vapour void fraction closure: dryness "
            "{:g}, mixture velocity {:g} m/s, drift flux velocity {:g} m/s, "
            "profile parameter {:g}.",
            dryness, v_mix, u_gu, C_0));
    }

    double const S = closureSlipParameter(state);

    // Upper bound of the admissible interval. The slip transports vapour out
    // of the control volume, so the void fraction stays below the homogeneous
    // one, alpha <= dryness / S <= alpha_homogeneous. That bound is at most
    // one for a profile parameter of at least one, which is what the
    // Rouhani-Axelsson correlation of driftFluxProfileParameter() gives. A
    // profile parameter below one describes a void profile peaking at the
    // wall rather than at the centre; the bound then exceeds one and the
    // volume fraction of the vapour, which cannot, hence the minimum.
    double const alpha_max = std::min(1., dryness / S);

    // The drift is aligned with the mixture flow, see above.
    auto const [a, b, c] = voidFractionQuadratic(alignDriftWithFlow(state));

    // The quadratic degenerates to a linear equation at vanishing mixture
    // velocity and at the critical point, where both phase densities coincide.
    // The tolerance is relative, to the linear coefficient below and to
    // b^2 for the discriminant, and some tens of machine epsilons, 45 of
    // them: that is the scale on which the coefficients, each a difference
    // of products of the closure state, lose their last digits, while a
    // quadratic whose leading coefficient is that much smaller than the linear
    // one has its small root within the accuracy of the linear solve anyway.
    constexpr double degeneracy_tolerance = 1e-14;

    // A root and the bound alpha_max are computed by different expressions, so
    // a root that coincides with the bound can come out just outside of it.
    // Roots are therefore accepted with a tolerance and clamped afterwards.
    // The tolerance is about the square root of the machine epsilon, which is
    // the accuracy a root of a quadratic is worth near a double root, where the
    // root shifts with the square root of the perturbation of the
    // coefficients.
    constexpr double interval_tolerance = 1e-8;

    auto const admissible = [&](double const alpha) -> std::optional<double>
    {
        double const tolerance = interval_tolerance * std::max(1., alpha_max);
        // The negated comparisons are deliberate; the plain form lets a NaN
        // root through to the clamp below, which returns it as an admissible
        // void fraction.
        if (!(alpha >= -tolerance) || !(alpha <= alpha_max + tolerance))
        {
            return std::nullopt;
        }
        return std::clamp(alpha, 0., alpha_max);
    };

    if (std::abs(a) <= degeneracy_tolerance * std::abs(b))
    {
        if (b == 0)
        {
            // Every coefficient vanishes, which for a validated two-phase
            // state happens exactly for a mixture at rest without drift, that
            // is v_mix = 0 and u_gu = 0. The closure then holds for every void
            // fraction. Retreating would not resolve it, because repeating the
            // time step does not change the velocity of the current iterate,
            // so the value is taken from the limit instead: without drift the
            // closure reads x = alpha S for every non-zero mixture velocity,
            // independently of it, hence alpha = x / S = alpha_max, the
            // profile slip alone.
            if (c == 0)
            {
                return admissible(alpha_max);
            }
            return std::nullopt;
        }
        return admissible(-c / b);
    }

    // Both terms underflow to zero for a mixture velocity in the subnormal
    // range, which loses the small root and is reported as no admissible root
    // rather than as a wrong one. Rescaling the coefficients would avoid it,
    // at the price of perturbing the rounding of every reachable state, and
    // the velocities in question are some 1e-160 m/s.
    double discriminant = b * b - 4 * a * c;
    if (discriminant < 0)
    {
        // A double root cannot be distinguished from a pair of complex roots
        // within the accuracy of the cancelling difference above.
        if (discriminant < -degeneracy_tolerance * b * b)
        {
            return std::nullopt;
        }
        discriminant = 0;
    }

    // Stable roots: the direct formula loses the small root to cancellation
    // whenever 4 a c is small compared to b^2, which is the case for small
    // dryness.
    double const q = -0.5 * (b + std::copysign(std::sqrt(discriminant), b));
    if (q == 0)
    {
        return admissible(0.);
    }

    // The smaller root is the one continuous with the single phase limit
    // alpha(dryness -> 0) = 0.
    double const root_1 = q / a;
    double const root_2 = c / q;

    auto const first = admissible(std::min(root_1, root_2));
    return first ? first : admissible(std::max(root_1, root_2));
}

double mixtureSlipParameter(double const alpha, DriftFluxState const& state)
{
    // The zero slip limit at a void fraction of one, see the documentation of
    // this function.
    if (alpha == 1)
    {
        return 0.;
    }

    double const rho_v = state.vapour_water_density;
    double const rho_l = state.liquid_water_density;
    double const C_0 = state.C_0;

    double const rho_mix = alpha * rho_v + (1 - alpha) * rho_l;

    return alpha * rho_l * rho_v * rho_mix / (1 - alpha) /
           std::pow((alpha * C_0 * rho_v + (1 - alpha * C_0) * rho_l), 2) *
           std::pow((C_0 - 1) * state.v_mix + state.u_gu, 2);
}
}  // namespace MaterialPropertyLib
