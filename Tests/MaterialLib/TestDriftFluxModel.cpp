// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>
#include <spdlog/fmt/fmt.h>

#include <cmath>
#include <limits>
#include <optional>
#include <string>

#include "MaterialLib/MPL/Utils/DriftFluxModel.h"
#include "MaterialLib/PhysicalConstant.h"
#include "NumLib/Exceptions.h"

namespace
{
// Saturated water and steam at about 1 MPa.
constexpr double rho_l = 887.1;
constexpr double rho_v = 5.145;

// The closure state at the saturated densities above. The tests that vary a
// density spell the state out instead.
MaterialPropertyLib::DriftFluxState saturatedState(double const dryness,
                                                   double const v_mix,
                                                   double const C_0,
                                                   double const u_gu)
{
    return {.dryness = dryness,
            .vapour_water_density = rho_v,
            .liquid_water_density = rho_l,
            .v_mix = v_mix,
            .C_0 = C_0,
            .u_gu = u_gu};
}

// Void fraction of a homogeneous mixture, i.e. without slip between the
// phases. It bounds the admissible interval from above, loosely: the bound the
// closure enforces is x / S, which is smaller whenever the profile parameter
// exceeds one, see admissibleUpperBound() below.
double homogeneousVoidFraction(double const dryness)
{
    return dryness * rho_l / (dryness * rho_l + (1 - dryness) * rho_v);
}

// Upper bound x / S of the admissible interval, cf. DriftFluxModel.h.
double admissibleUpperBound(double const dryness)
{
    double const C_0 = MaterialPropertyLib::driftFluxProfileParameter(dryness);
    return dryness / (C_0 * (dryness + (1 - dryness) * rho_v / rho_l));
}
}  // namespace

// The quadratic coefficients must reproduce the drift-flux residual for
// arbitrary void fractions, not only in the roots.
TEST(MaterialLibDriftFluxModel, QuadraticMatchesResidual)
{
    using namespace MaterialPropertyLib;

    double const dryness = 0.3;
    double const C_0 = MaterialPropertyLib::driftFluxProfileParameter(dryness);
    double const u_gu = 0.21;
    double const v_mix = 1.7;

    auto const state = saturatedState(dryness, v_mix, C_0, u_gu);
    auto const q = voidFractionQuadratic(state);

    for (double alpha = -0.5; alpha <= 1.5; alpha += 0.1)
    {
        double const from_quadratic =
            rho_l * ((q.a * alpha + q.b) * alpha + q.c);
        double const residual = voidFractionResidual(alpha, state);
        EXPECT_NEAR(residual, from_quadratic, 1e-6 * std::abs(rho_l * q.c))
            << "alpha = " << alpha;
    }
}

// For co-current upflow the residual changes sign on (0, 1), hence a root
// exists for every dryness and every positive mixture velocity.
TEST(MaterialLibDriftFluxModel, UpflowRootExistsAndIsAdmissible)
{
    using namespace MaterialPropertyLib;

    double const u_gu = 0.21;

    for (double dryness = 0.01; dryness < 1.; dryness += 0.01)
    {
        double const C_0 =
            MaterialPropertyLib::driftFluxProfileParameter(dryness);

        for (double v_mix : {1e-4, 1e-2, 0.5, 1.7, 12., 150.})
        {
            auto const state = saturatedState(dryness, v_mix, C_0, u_gu);
            auto const alpha = computeVapourVoidFraction(state);

            ASSERT_TRUE(alpha.has_value())
                << "dryness = " << dryness << ", v = " << v_mix;
            EXPECT_GE(*alpha, 0.) << "dryness = " << dryness;
            EXPECT_LE(*alpha, admissibleUpperBound(dryness))
                << "dryness = " << dryness << ", v = " << v_mix;
            EXPECT_LE(*alpha, homogeneousVoidFraction(dryness));

            double const residual = voidFractionResidual(*alpha, state);
            double const scale = std::abs(voidFractionResidual(0., state));
            EXPECT_NEAR(residual / scale, 0., 1e-10)
                << "dryness = " << dryness << ", v = " << v_mix;
        }
    }
}

// Without slip (no drift velocity, unit profile parameter) the closure must
// return the homogeneous void fraction.
TEST(MaterialLibDriftFluxModel, NoSlipGivesHomogeneousMixture)
{
    using namespace MaterialPropertyLib;

    for (double dryness = 0.05; dryness < 1.; dryness += 0.05)
    {
        auto const alpha =
            computeVapourVoidFraction(saturatedState(dryness, 1.7, 1., 0.));

        ASSERT_TRUE(alpha.has_value());
        EXPECT_NEAR(*alpha, homogeneousVoidFraction(dryness), 1e-12);
    }
}

// Single phase states are handled without solving the quadratic.
TEST(MaterialLibDriftFluxModel, SinglePhaseEndPoints)
{
    using namespace MaterialPropertyLib;

    EXPECT_EQ(0.,
              computeVapourVoidFraction(saturatedState(0., 1.7, 1.12, 0.21)));
    EXPECT_EQ(1., computeVapourVoidFraction(saturatedState(1., 1.7, 1., 0.)));
}

// A non-finite state follows a diverging global Newton step, so the assembly
// is aborted instead of the run being ended. The Newton solver turns that into
// a repeated time step.
TEST(MaterialLibDriftFluxModel, NonFiniteStateAbortsAssembly)
{
    using namespace MaterialPropertyLib;

    double const nan = std::numeric_limits<double>::quiet_NaN();
    double const inf = std::numeric_limits<double>::infinity();

    EXPECT_THROW(
        computeVapourVoidFraction(saturatedState(0.3, nan, 1.12, 0.21)),
        NumLib::AssemblyException);
    EXPECT_THROW(
        computeVapourVoidFraction(saturatedState(0.3, inf, 1.12, 0.21)),
        NumLib::AssemblyException);
    EXPECT_THROW(computeVapourVoidFraction(saturatedState(0.3, 1.7, 1.12, nan)),
                 NumLib::AssemblyException);
    EXPECT_THROW(computeVapourVoidFraction(saturatedState(0.3, 1.7, nan, 0.21)),
                 NumLib::AssemblyException);
}

// The saturation line correlations are extrapolated outside of their pressure
// range and then return non-positive densities, which aborts the assembly.
TEST(MaterialLibDriftFluxModel, NonPositiveDensityAbortsAssembly)
{
    using namespace MaterialPropertyLib;

    for (double vapour_density : {-25.8, 0.})
    {
        EXPECT_THROW(
            computeVapourVoidFraction({.dryness = 0.3,
                                       .vapour_water_density = vapour_density,
                                       .liquid_water_density = rho_l,
                                       .v_mix = 1.7,
                                       .C_0 = 1.12,
                                       .u_gu = 0.21}),
            NumLib::AssemblyException)
            << "vapour density = " << vapour_density;
    }
    for (double liquid_density : {-502.5, 0.})
    {
        EXPECT_THROW(
            computeVapourVoidFraction({.dryness = 0.3,
                                       .vapour_water_density = rho_v,
                                       .liquid_water_density = liquid_density,
                                       .v_mix = 1.7,
                                       .C_0 = 1.12,
                                       .u_gu = 0.21}),
            NumLib::AssemblyException)
            << "liquid density = " << liquid_density;
    }
}

// A NaN density is what the extrapolated correlations return where the base of
// a non-integer power turns negative, and it is the reason the guard is
// written as a negated comparison: `density <= 0` would let it through into
// the closure, which would then hand out a NaN void fraction.
TEST(MaterialLibDriftFluxModel, NonFiniteDensityAbortsAssembly)
{
    using namespace MaterialPropertyLib;

    double const nan = std::numeric_limits<double>::quiet_NaN();
    double const inf = std::numeric_limits<double>::infinity();

    for (double vapour_density : {nan, -inf})
    {
        EXPECT_THROW(
            computeVapourVoidFraction({.dryness = 0.3,
                                       .vapour_water_density = vapour_density,
                                       .liquid_water_density = rho_l,
                                       .v_mix = 1.7,
                                       .C_0 = 1.12,
                                       .u_gu = 0.21}),
            NumLib::AssemblyException)
            << "vapour density = " << vapour_density;
    }
    EXPECT_THROW(computeVapourVoidFraction({.dryness = 0.3,
                                            .vapour_water_density = rho_v,
                                            .liquid_water_density = nan,
                                            .v_mix = 1.7,
                                            .C_0 = 1.12,
                                            .u_gu = 0.21}),
                 NumLib::AssemblyException);
}

// In a single phase state the density of the absent phase is the one
// extrapolated far off the saturation line, and it does not enter the result,
// so it must not abort the assembly.
TEST(MaterialLibDriftFluxModel, SinglePhaseIgnoresAbsentPhaseDensity)
{
    using namespace MaterialPropertyLib;

    double const nan = std::numeric_limits<double>::quiet_NaN();

    for (double vapour_density : {-25.8, nan})
    {
        EXPECT_EQ(0., computeVapourVoidFraction(
                          {.dryness = 0.,
                           .vapour_water_density = vapour_density,
                           .liquid_water_density = rho_l,
                           .v_mix = 1.7,
                           .C_0 = 1.12,
                           .u_gu = 0.21}));
    }
    for (double liquid_density : {-502.5, nan})
    {
        EXPECT_EQ(1., computeVapourVoidFraction(
                          {.dryness = 1.,
                           .vapour_water_density = rho_v,
                           .liquid_water_density = liquid_density,
                           .v_mix = 1.7,
                           .C_0 = 1.,
                           .u_gu = 0.}));
    }
}

// A NaN dryness passes neither single phase comparison and must not reach the
// closure, where it would come back as a NaN void fraction.
TEST(MaterialLibDriftFluxModel, NonFiniteDrynessAbortsAssembly)
{
    using namespace MaterialPropertyLib;

    double const nan = std::numeric_limits<double>::quiet_NaN();

    EXPECT_THROW(
        computeVapourVoidFraction(saturatedState(nan, 1.7, 1.12, 0.21)),
        NumLib::AssemblyException);
}

// Slip reduces the void fraction monotonically: a larger drift velocity lets
// more vapour leave the control volume, so less of it is present.
TEST(MaterialLibDriftFluxModel, DriftReducesVoidFraction)
{
    using namespace MaterialPropertyLib;

    double const dryness = 0.3;
    double const C_0 = MaterialPropertyLib::driftFluxProfileParameter(dryness);

    double previous = homogeneousVoidFraction(dryness);
    for (double u_gu : {0., 0.05, 0.21, 0.8, 3.})
    {
        auto const alpha =
            computeVapourVoidFraction(saturatedState(dryness, 1.7, C_0, u_gu));
        ASSERT_TRUE(alpha.has_value());
        EXPECT_LE(*alpha, previous) << "u_gu = " << u_gu;
        previous = *alpha;
    }
}

// The quadratic degenerates to a linear equation when the phase densities
// coincide (critical point) or the mixture is at rest.
TEST(MaterialLibDriftFluxModel, DegenerateQuadratic)
{
    using namespace MaterialPropertyLib;

    double const dryness = 0.3;
    double const C_0 = MaterialPropertyLib::driftFluxProfileParameter(dryness);

    // Equal densities: void fraction equals dryness up to the slip.
    {
        DriftFluxState const state{.dryness = dryness,
                                   .vapour_water_density = rho_l,
                                   .liquid_water_density = rho_l,
                                   .v_mix = 1.7,
                                   .C_0 = C_0,
                                   .u_gu = 0.21};

        auto const alpha = computeVapourVoidFraction(state);
        ASSERT_TRUE(alpha.has_value());
        EXPECT_GE(*alpha, 0.);
        EXPECT_LE(*alpha, dryness / C_0);
        EXPECT_NEAR(voidFractionResidual(*alpha, state), 0., 1e-6);
    }

    // Mixture at rest: all vapour drifts away, no vapour is retained.
    {
        auto const alpha =
            computeVapourVoidFraction(saturatedState(dryness, 0., C_0, 0.21));
        ASSERT_TRUE(alpha.has_value());
        EXPECT_NEAR(*alpha, 0., 1e-12);
    }
}

// Where the two phase densities cross, in the immediate vicinity of the
// critical point, the parabola opens the other way and the smaller of its two
// roots leaves the admissible interval, so the larger one is the physical one.
// This is the reason the root is selected by testing the interval instead of
// by its position, and it is the only state that exercises that selection.
TEST(MaterialLibDriftFluxModel, CrossedDensitiesTakeTheSecondRoot)
{
    using namespace MaterialPropertyLib;

    double const dryness = 0.1;
    double const C_0 = driftFluxProfileParameter(dryness);

    // Beyond the critical point the vapour is the denser phase.
    double const vapour_density = 1.0001 * rho_l;
    double const bound =
        dryness / (C_0 * (dryness + (1 - dryness) * vapour_density / rho_l));

    for (double const v_mix : {-1e3, -1e5})
    {
        DriftFluxState const state{.dryness = dryness,
                                   .vapour_water_density = vapour_density,
                                   .liquid_water_density = rho_l,
                                   .v_mix = v_mix,
                                   .C_0 = C_0,
                                   .u_gu = 0.21};

        auto const alpha = computeVapourVoidFraction(state);

        ASSERT_TRUE(alpha.has_value()) << "v = " << v_mix;
        EXPECT_GE(*alpha, 0.) << "v = " << v_mix;
        EXPECT_LE(*alpha, bound) << "v = " << v_mix;

        // The other root of the same quadratic is far outside the interval,
        // which is what makes this state the one that selects the second. The
        // quadratic is the one of the aligned state, which is what
        // computeVapourVoidFraction() solves.
        DriftFluxState aligned = state;
        aligned.u_gu = std::copysign(state.u_gu, state.v_mix);

        auto const q = voidFractionQuadratic(aligned);
        double const other_root = q.c / (q.a * *alpha);
        EXPECT_LT(other_root, -1.) << "v = " << v_mix;

        // The crossed state is ill-conditioned: the density difference is
        // four orders of magnitude below the densities themselves, so the
        // quadratic is nearly degenerate and its root carries a correspondingly
        // larger residual than the well separated states above.
        double const residual = voidFractionResidual(*alpha, state);
        double const scale = std::abs(voidFractionResidual(0., state));
        EXPECT_NEAR(residual / scale, 0., 1e-3) << "v = " << v_mix;
    }
}

// Existence and uniqueness hold for the closure at every mixture velocity, but
// not for the arithmetic that solves it: below some 1e-160 m/s the
// coefficients of the quadratic underflow and the small root is lost. The
// documented answer is the same as for any state without an admissible root,
// no value, which the assembly turns into a retreat.
TEST(MaterialLibDriftFluxModel, SubnormalMixtureVelocityHasNoAdmissibleRoot)
{
    using namespace MaterialPropertyLib;

    double const dryness = 0.3;
    double const C_0 = MaterialPropertyLib::driftFluxProfileParameter(dryness);

    // Still solvable an order of magnitude above the underflow.
    EXPECT_TRUE(
        computeVapourVoidFraction(saturatedState(dryness, 1e-150, C_0, 0.))
            .has_value());

    for (double const v_mix : {1e-200, -1e-200, 1e-300})
    {
        EXPECT_FALSE(
            computeVapourVoidFraction(saturatedState(dryness, v_mix, C_0, 0.))
                .has_value())
            << "v = " << v_mix;
    }
}

// A mixture at rest without drift satisfies the closure with every void
// fraction: every coefficient of the quadratic vanishes. Without drift the
// closure is solved by the profile slip alone for every non-zero mixture
// velocity, independently of it, so that value is also the one at rest.
TEST(MaterialLibDriftFluxModel, RestWithoutDriftIsTheProfileSlipLimit)
{
    using namespace MaterialPropertyLib;

    double const dryness = 0.3;
    double const C_0 = driftFluxProfileParameter(dryness);

    auto const state = saturatedState(dryness, 0., C_0, 0.);

    // The residual vanishes for every void fraction, so nothing in the state
    // at rest singles one out.
    EXPECT_EQ(0., voidFractionResidual(0., state));
    EXPECT_EQ(0., voidFractionResidual(0.5, state));
    EXPECT_EQ(0., voidFractionResidual(1., state));

    auto const at_rest = computeVapourVoidFraction(state);
    ASSERT_TRUE(at_rest.has_value());
    EXPECT_DOUBLE_EQ(admissibleUpperBound(dryness), *at_rest);

    // It is the limit of the drift free closure, which is independent of the
    // mixture velocity, and it is approached from both flow directions. At
    // rest the value is x / S itself, whereas for a moving mixture it is a
    // root of the quadratic, so the two agree only up to round-off.
    for (double v_mix : {-25., -1.33, -1e-6, 1e-6, 1.33, 25.})
    {
        auto const moving =
            computeVapourVoidFraction(saturatedState(dryness, v_mix, C_0, 0.));
        ASSERT_TRUE(moving.has_value()) << "v = " << v_mix;
        EXPECT_NEAR(1., *moving / *at_rest, 1e-13) << "v = " << v_mix;
    }
}

// A vanishing mixture velocity with a drift is not degenerate: the drift
// carries all vapour out of the control volume, so the void fraction is zero.
TEST(MaterialLibDriftFluxModel, RestWithDriftHasNoVapour)
{
    using namespace MaterialPropertyLib;

    double const dryness = 0.3;
    double const C_0 = driftFluxProfileParameter(dryness);

    EXPECT_EQ(
        0., computeVapourVoidFraction(saturatedState(dryness, 0., C_0, 0.21)));
}

// Small dryness, large density ratio: the stable quadratic formula must not
// lose the small root to cancellation.
TEST(MaterialLibDriftFluxModel, SmallDrynessIsAccurate)
{
    using namespace MaterialPropertyLib;

    double const u_gu = 0.21;
    double const v_mix = 1.7;

    for (double dryness : {1e-3, 1e-6, 1e-9, 1e-12})
    {
        double const C_0 =
            MaterialPropertyLib::driftFluxProfileParameter(dryness);
        auto const alpha = computeVapourVoidFraction(
            saturatedState(dryness, v_mix, C_0, u_gu));

        ASSERT_TRUE(alpha.has_value()) << "dryness = " << dryness;
        EXPECT_GT(*alpha, 0.) << "dryness = " << dryness;
        EXPECT_LE(*alpha, admissibleUpperBound(dryness))
            << "dryness = " << dryness;

        // The residual is linear in the dryness for vanishing dryness, so the
        // void fraction is too. Its slope is the small root of the linearised
        // closure, alpha / x -> v rho_l / (rho_v (C_0 v + u_gu)), which stays
        // below the homogeneous slope 1 / S because the drift carries vapour
        // away. The relative deviation from the slope is of the order of
        // 170 x, so the tolerance is taken proportional to the dryness; a
        // fixed band would not distinguish the two slopes, which differ by
        // ten per cent.
        double const slope =
            v_mix * rho_l /
            (rho_v *
             (MaterialPropertyLib::driftFluxProfileParameter(0.) * v_mix +
              u_gu));
        EXPECT_NEAR(*alpha / dryness, slope, slope * (200 * dryness + 1e-12))
            << "dryness = " << dryness;
    }
}

// Backflow: the drift is aligned with the mixture flow, so the void fraction
// is an even function of the mixture velocity and stays admissible.
TEST(MaterialLibDriftFluxModel, BackflowMirrorsUpflow)
{
    using namespace MaterialPropertyLib;

    double const u_gu = 0.21;

    for (double dryness = 0.01; dryness < 1.; dryness += 0.01)
    {
        double const C_0 =
            MaterialPropertyLib::driftFluxProfileParameter(dryness);

        for (double v_mix : {1e-3, 0.1, 1.33, 25.})
        {
            auto const upflow = computeVapourVoidFraction(
                saturatedState(dryness, v_mix, C_0, u_gu));
            auto const backflow = computeVapourVoidFraction(
                saturatedState(dryness, -v_mix, C_0, u_gu));

            ASSERT_TRUE(backflow.has_value())
                << "dryness = " << dryness << ", v = " << -v_mix;
            ASSERT_TRUE(upflow.has_value());
            EXPECT_DOUBLE_EQ(*upflow, *backflow)
                << "dryness = " << dryness << ", v = " << v_mix;
            EXPECT_LE(*backflow, admissibleUpperBound(dryness));
        }
    }
}

// The alignment only sets the sign, so it leaves the magnitude of the drift
// flux velocity alone and passing an already aligned value through it again
// must not change it.
TEST(MaterialLibDriftFluxModel, DriftAlignmentIsIdempotent)
{
    using namespace MaterialPropertyLib;

    double const u_gu = 0.21;

    for (double v_mix : {-25., -1.33, -1e-3, 0., 1e-3, 1.33, 25.})
    {
        double const aligned = alignedDriftFluxVelocity(u_gu, v_mix);

        EXPECT_EQ(u_gu, std::abs(aligned)) << "v = " << v_mix;
        EXPECT_EQ(aligned, alignedDriftFluxVelocity(aligned, v_mix))
            << "v = " << v_mix;
    }
}

// The slip momentum term of the mixture is built from the same bracket
// (C_0 - 1) v + u_gu that the closure solves with, so with the aligned drift
// flux velocity it is mirror symmetric in the mixture velocity, just as the
// void fraction is. With a raw drift flux velocity it would instead vanish at
// the finite backflow velocity v = -u_gu / (C_0 - 1).
TEST(MaterialLibDriftFluxModel, SlipBracketMirrorsUnderFlowReversal)
{
    using namespace MaterialPropertyLib;

    double const u_gu = 0.21;
    double const dryness = 0.3;
    double const C_0 = driftFluxProfileParameter(dryness);

    auto const slipBracket = [&](double const v)
    { return (C_0 - 1) * v + alignedDriftFluxVelocity(u_gu, v); };

    for (double v_mix : {1e-3, 0.1, 1.33, 25.})
    {
        EXPECT_DOUBLE_EQ(slipBracket(v_mix), -slipBracket(-v_mix))
            << "v = " << v_mix;
    }

    EXPECT_NE(0., slipBracket(-u_gu / (C_0 - 1)));
}

// The void fraction is continuous at flow reversal, where it vanishes
// proportionally to the magnitude of the mixture velocity.
TEST(MaterialLibDriftFluxModel, ContinuousAtFlowReversal)
{
    using namespace MaterialPropertyLib;

    double const dryness = 0.3;
    double const C_0 = MaterialPropertyLib::driftFluxProfileParameter(dryness);
    double const u_gu = 0.21;

    double previous = std::numeric_limits<double>::max();
    for (double v_mix : {1e-2, 1e-4, 1e-6, 1e-8, 1e-10})
    {
        auto const alpha = computeVapourVoidFraction(
            saturatedState(dryness, v_mix, C_0, u_gu));
        ASSERT_TRUE(alpha.has_value()) << "v = " << v_mix;
        EXPECT_GT(*alpha, 0.) << "v = " << v_mix;
        EXPECT_LT(*alpha, previous) << "v = " << v_mix;
        previous = *alpha;
    }

    EXPECT_EQ(
        0., computeVapourVoidFraction(saturatedState(dryness, 0., C_0, u_gu)));
}

// Beyond the critical point the surface tension correlation and the buoyancy
// have negative bases, whose non-integer powers would be NaN. Both are capped,
// so the drift vanishes there instead.
TEST(MaterialLibDriftFluxModel, DriftVanishesBeyondCriticalPoint)
{
    using namespace MaterialPropertyLib;

    double const T_critical =
        MaterialLib::PhysicalConstant::CriticalPoint::TemperatureWater;

    EXPECT_GT(driftFluxVelocity(0.3, 453.03, rho_v, rho_l), 0.);
    EXPECT_EQ(0., driftFluxVelocity(0.3, T_critical, rho_v, rho_l));
    EXPECT_EQ(0., driftFluxVelocity(0.3, T_critical + 50., rho_v, rho_l));

    // Below the critical point the densities can still cross, which the
    // buoyancy cap catches separately from the surface tension.
    EXPECT_EQ(0., driftFluxVelocity(0.3, 453.03, rho_l, rho_v));
}

// The caps at the critical point must not swallow a non-finite temperature: a
// drift flux velocity of zero would look like an admissible state to the
// closure, whereas the NaN aborts the assembly.
TEST(MaterialLibDriftFluxModel, NonFiniteTemperatureAbortsAssembly)
{
    using namespace MaterialPropertyLib;

    double const nan = std::numeric_limits<double>::quiet_NaN();
    double const dryness = 0.3;
    double const u_gu = driftFluxVelocity(dryness, nan, rho_v, rho_l);

    EXPECT_TRUE(std::isnan(u_gu));
    EXPECT_THROW(
        computeVapourVoidFraction(saturatedState(
            dryness, 1.7,
            MaterialPropertyLib::driftFluxProfileParameter(dryness), u_gu)),
        NumLib::AssemblyException);
}

// The drift vanishes in the single phase vapour limit, which is what makes the
// slip parameter of the mixture regular at a void fraction of one.
TEST(MaterialLibDriftFluxModel, DriftVanishesAtFullDryness)
{
    using namespace MaterialPropertyLib;

    EXPECT_EQ(0., driftFluxVelocity(1., 453.03, rho_v, rho_l));
    EXPECT_GT(driftFluxVelocity(0.99, 453.03, rho_v, rho_l), 0.);
}

// The two powers of the drift flux correlation, the fourth root of the
// buoyancy and the inverse square root of the liquid density, are what set the
// magnitude of the slip, so they are pinned by their scaling rather than by a
// reproduced number: sixteen times the buoyancy doubles the drift, and four
// times the liquid density at unchanged buoyancy halves it.
TEST(MaterialLibDriftFluxModel, DriftScalesWithTheCorrelationPowers)
{
    using namespace MaterialPropertyLib;

    double const dryness = 0.3;
    double const temperature = 453.03;

    // Buoyancy is proportional to the density difference, so a difference of
    // rho_l against one of rho_l / 16 is a factor of sixteen under the fourth
    // root, that is a factor of two.
    double const drift_full =
        driftFluxVelocity(dryness, temperature, 0., rho_l);
    double const drift_sixteenth =
        driftFluxVelocity(dryness, temperature, rho_l - rho_l / 16., rho_l);

    EXPECT_NEAR(2., drift_full / drift_sixteenth, 1e-12);

    // The liquid density enters the buoyancy as well, so it is held at the
    // same difference while the density itself is quadrupled: only the inverse
    // square root remains, that is a factor of one half.
    double const difference = rho_l - rho_v;
    double const drift_dense = driftFluxVelocity(
        dryness, temperature, 4 * rho_l - difference, 4 * rho_l);

    EXPECT_NEAR(
        0.5,
        drift_dense / driftFluxVelocity(dryness, temperature, rho_v, rho_l),
        1e-12);
}

// The scaling test above holds the temperature fixed, so the surface tension
// factor of the correlation is the same on both sides of every ratio it forms
// and cancels out of all of them: its three constants could take any value
// whatsoever and it would still pass. They are the IAPWS R1-76 correlation,
// sigma = 235.8 mN/m * tau^1.256 * (1 - 0.625 tau) with tau = 1 - T / T_c, so
// they are pinned here against that correlation's own published values rather
// than against a number reproduced from this implementation.
TEST(MaterialLibDriftFluxModel, DriftCarriesTheSurfaceTensionCorrelation)
{
    using namespace MaterialPropertyLib;

    // At vanishing dryness and vapour density the drift velocity reduces to
    // 1.18 (g sigma rho_l)^(1/4) / sqrt(rho_l), which inverts to the surface
    // tension that went into it.
    auto const surface_tension = [](double const temperature)
    {
        constexpr double gravity = 9.81;          // m/s^2, as in the closure
        constexpr double density_liquid = 1000.;  // kg/m^3

        double const drift =
            driftFluxVelocity(0., temperature, 0., density_liquid);

        return std::pow(drift * std::sqrt(density_liquid) / 1.18, 4) /
               (gravity * density_liquid);
    };

    // IAPWS R1-76: 71.69 mN/m at 300 K and 58.91 mN/m at 373.15 K. The
    // tolerance is a hundredth of a millinewton per metre, which is the
    // precision the released values are quoted to; the exponent alone moves
    // the first of them by ten millinewtons per metre if it is wrong.
    EXPECT_NEAR(71.686e-3, surface_tension(300.), 1e-5);
    EXPECT_NEAR(58.912e-3, surface_tension(373.15), 1e-5);
}

// Without drift the closure reads x = alpha S, so the void fraction is the
// profile slip alone and equals the upper bound of the admissible interval. It
// must not come out above that bound: the root and the bound are computed by
// different expressions, so the root can exceed it by an ulp, and a void
// fraction above x / S is not a state the closure admits.
TEST(MaterialLibDriftFluxModel, WithoutDriftTheBoundIsAttainedNotExceeded)
{
    using namespace MaterialPropertyLib;

    for (double const dryness : {1e-9, 1e-3, 0.3, 0.9})
    {
        double const C_0 = driftFluxProfileParameter(dryness);

        for (double const density_ratio : {1e-4, 1e-3, 0.1})
        {
            double const vapour_density = density_ratio * rho_l;

            // The bound of the admissible interval, x / S, for these
            // densities; admissibleUpperBound() above is tied to the
            // saturated pair of the other tests.
            double const bound =
                dryness /
                (C_0 * (dryness + (1 - dryness) * vapour_density / rho_l));

            for (double const v_mix : {-100., -1.7, 1.7, 100.})
            {
                auto const alpha = computeVapourVoidFraction(
                    {.dryness = dryness,
                     .vapour_water_density = vapour_density,
                     .liquid_water_density = rho_l,
                     .v_mix = v_mix,
                     .C_0 = C_0,
                     .u_gu = 0.});

                ASSERT_TRUE(alpha.has_value())
                    << "x = " << dryness
                    << ", rho_v / rho_l = " << density_ratio
                    << ", v = " << v_mix;
                EXPECT_LE(*alpha, bound)
                    << "x = " << dryness
                    << ", rho_v / rho_l = " << density_ratio
                    << ", v = " << v_mix;
                EXPECT_NEAR(*alpha, bound, 1e-12 * bound)
                    << "x = " << dryness
                    << ", rho_v / rho_l = " << density_ratio
                    << ", v = " << v_mix;
            }
        }
    }
}

// The diagnostics function's contract is its string, so it is pinned here:
// the residuals it reports must be the ones at the ends of the interval that
// computeVapourVoidFraction() solves on, [0, x/S], and not at a void fraction
// of one. The second root of the closure is always larger than one, so a sign
// change over [0, 1] does not imply an admissible root, and reporting the
// residual there overstates what the numbers show.
TEST(MaterialLibDriftFluxModel, DiagnosticsReportTheAdmissibleInterval)
{
    using namespace MaterialPropertyLib;

    double const dryness = 0.3;
    double const C_0 = driftFluxProfileParameter(dryness);
    double const u_gu = 0.21;
    double const alpha_max = admissibleUpperBound(dryness);

    for (double const v_mix : {1.7, -1.7})
    {
        auto const message = voidFractionClosureDiagnostics(
            saturatedState(dryness, v_mix, C_0, u_gu));

        // The drift is aligned with the mixture flow, so the residuals of the
        // message are the ones of the closure that is actually solved.
        auto const aligned = saturatedState(
            dryness, v_mix, C_0, alignedDriftFluxVelocity(u_gu, v_mix));
        auto const residual_at = [&](double const alpha)
        { return fmt::format("{:g}", voidFractionResidual(alpha, aligned)); };

        EXPECT_NE(std::string::npos,
                  message.find(fmt::format("{:g}", alpha_max)))
            << "v = " << v_mix << ": " << message;
        EXPECT_NE(std::string::npos, message.find(residual_at(0.)))
            << "v = " << v_mix << ": " << message;
        EXPECT_NE(std::string::npos, message.find(residual_at(alpha_max)))
            << "v = " << v_mix << ": " << message;
        EXPECT_EQ(std::string::npos, message.find(residual_at(1.)))
            << "v = " << v_mix << ": " << message;
    }
}

// A profile parameter below one, which the Rouhani-Axelsson correlation does
// not produce but the closure accepts, lifts the bound x / S of the admissible
// interval above one. The void fraction is a volume fraction and must stay
// below one anyway, otherwise the liquid fraction 1 - alpha turns negative and
// the slip parameter of the mixture comes back negative with it.
TEST(MaterialLibDriftFluxModel, ProfileParameterBelowOneStaysBelowFullVoid)
{
    using namespace MaterialPropertyLib;

    double const v_mix = 1.7;
    double const u_gu = 0.21;

    for (double C_0 : {0.5, 0.8, 0.95})
    {
        for (double dryness : {0.3, 0.6, 0.9, 0.99})
        {
            auto const state = saturatedState(dryness, v_mix, C_0, u_gu);
            auto const alpha = computeVapourVoidFraction(state);

            // Whether such a state has an admissible root at all depends on
            // how far above one the bound x / S is lifted. What must hold
            // either way is that no void fraction above one is returned.
            if (alpha)
            {
                EXPECT_GE(*alpha, 0.)
                    << "C_0 = " << C_0 << ", dryness = " << dryness;
                EXPECT_LE(*alpha, 1.)
                    << "C_0 = " << C_0 << ", dryness = " << dryness;
                EXPECT_GE(mixtureSlipParameter(*alpha, state), 0.)
                    << "C_0 = " << C_0 << ", dryness = " << dryness;
            }
        }
    }

    // A state whose only root lies above one is refused rather than clamped to
    // one: the closure does not predict a pure vapour section there, and
    // saying so would hide an unusable profile parameter behind a plausible
    // looking void fraction. Without the bound this state returns 1.000954.
    EXPECT_FALSE(
        computeVapourVoidFraction(saturatedState(0.3, v_mix, 0.5, u_gu))
            .has_value());
}

// The slip parameter of the mixture vanishes in both single phase limits: at a
// void fraction of zero there is no vapour to slip, and at a void fraction of
// one no liquid to slip against, where the expression would divide by zero.
TEST(MaterialLibDriftFluxModel, MixtureSlipVanishesInSinglePhaseLimits)
{
    using namespace MaterialPropertyLib;

    double const v_mix = 1.7;
    double const u_gu = 0.21;

    EXPECT_EQ(
        0.,
        mixtureSlipParameter(
            0., saturatedState(
                    0., v_mix,
                    MaterialPropertyLib::driftFluxProfileParameter(0.), u_gu)));
    EXPECT_EQ(
        0.,
        mixtureSlipParameter(
            1., saturatedState(
                    1., v_mix,
                    MaterialPropertyLib::driftFluxProfileParameter(1.), 0.)));
    // The limit holds whatever the rest of the state is: at a void fraction of
    // one the expression is not evaluated at all.
    EXPECT_EQ(0.,
              mixtureSlipParameter(
                  1., saturatedState(
                          0.3, v_mix,
                          MaterialPropertyLib::driftFluxProfileParameter(0.3),
                          u_gu)));
}

// The slip parameter approaches its zero limit continuously from below a void
// fraction of one, along the physical path where the drift vanishes with the
// liquid fraction, and it does so linearly in 1 - dryness. That order is what
// makes the limit finite at all: the bracket (C_0 - 1) v + u_gu enters
// squared and the liquid fraction 1 - alpha divides, and the two cancel to
// first order. A quadratic or constant decay would mean one of the two factors
// does not vanish as documented.
TEST(MaterialLibDriftFluxModel, MixtureSlipVanishesLinearlyAtFullVoidFraction)
{
    using namespace MaterialPropertyLib;

    double const v_mix = 1.7;
    double const temperature = 453.03;

    double previous_gamma = std::numeric_limits<double>::max();
    double previous_slope = 0.;
    for (double liquid_fraction : {1e-2, 1e-3, 1e-4, 1e-5})
    {
        double const dryness = 1 - liquid_fraction;
        double const C_0 =
            MaterialPropertyLib::driftFluxProfileParameter(dryness);
        double const u_gu =
            driftFluxVelocity(dryness, temperature, rho_v, rho_l);

        auto const state = saturatedState(dryness, v_mix, C_0, u_gu);
        auto const alpha = computeVapourVoidFraction(state);
        ASSERT_TRUE(alpha.has_value()) << "dryness = " << dryness;

        double const gamma = mixtureSlipParameter(*alpha, state);
        EXPECT_GT(gamma, 0.) << "dryness = " << dryness;
        EXPECT_LT(gamma, previous_gamma) << "dryness = " << dryness;
        previous_gamma = gamma;

        // The slope of the linear decay, which settles at about 554 for this
        // state. Only the two finest steps are compared, where the higher
        // order terms have died out; the coarse ones still carry them.
        double const slope = gamma / liquid_fraction;
        if (liquid_fraction <= 1e-5)
        {
            EXPECT_NEAR(slope, previous_slope, 0.01 * slope)
                << "dryness = " << dryness;
        }
        previous_slope = slope;
    }
}

// The state factory is the single place the closure's inputs are composed, so
// what it returns has to be the state the local assemblers built by hand
// before: the profile parameter of the dryness, and the drift flux velocity
// aligned with the mixture flow rather than the raw one.
TEST(MaterialLibDriftFluxModel, StateFactoryComposesTheClosureInputs)
{
    using namespace MaterialPropertyLib;

    double const dryness = 0.3;
    double const temperature = 453.03;

    for (double v_mix : {-25., -1.33, 0., 1.33, 25.})
    {
        DriftFluxState const state =
            driftFluxState(dryness, temperature, rho_v, rho_l, v_mix);

        EXPECT_EQ(dryness, state.dryness) << "v = " << v_mix;
        EXPECT_EQ(rho_v, state.vapour_water_density) << "v = " << v_mix;
        EXPECT_EQ(rho_l, state.liquid_water_density) << "v = " << v_mix;
        EXPECT_EQ(v_mix, state.v_mix) << "v = " << v_mix;
        EXPECT_EQ(driftFluxProfileParameter(dryness), state.C_0)
            << "v = " << v_mix;
        EXPECT_EQ(
            alignedDriftFluxVelocity(
                driftFluxVelocity(dryness, temperature, rho_v, rho_l), v_mix),
            state.u_gu)
            << "v = " << v_mix;
    }
}
