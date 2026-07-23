// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <string>
#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"
#include "NumLib/ODESolver/AndersonAcceleration.h"
#include "NumLib/ODESolver/AndersonWeights.h"

// Note: no file-scope using-declaration - test sources are compiled in unity
// builds, where it would leak into the whole blob.

namespace
{
//! A Gram matrix together with the weights it must produce.
struct AndersonWeightsCase
{
    std::string name;
    Eigen::MatrixXd gram;
    Eigen::VectorXd expected_theta;
    double tolerance;
};

//! Builds the Gram matrix \f$ G = F^T F \f$ of the residuals given as the
//! columns of \p F, i.e. the input NumLib::detail::computeAndersonWeights() is
//! documented for.
Eigen::MatrixXd gramOf(Eigen::MatrixXd const& F)
{
    return F.transpose() * F;
}

std::vector<AndersonWeightsCase> exactWeightsCases()
{
    std::vector<AndersonWeightsCase> cases;

    {
        // Two residuals of equal norm that are orthogonal to each other:
        // G = ||f||^2 * I. By symmetry the norm-minimizing mixture under the
        // sum-to-one constraint must weight both equally.
        Eigen::MatrixXd G(2, 2);
        G << 2.0, 0.0,  //
            0.0, 2.0;
        cases.push_back({"TwoEqualOrthogonalResidualsSplitEvenly", G,
                         Eigen::Vector2d(0.5, 0.5), 1e-12});
    }
    {
        // f_0 is the zero vector (already converged), f_1 is not. Since
        // ||sum theta_i f_i|| is minimized and only f_0 can contribute zero
        // norm while satisfying sum(theta) = 1, the optimum is theta = [1, 0].
        Eigen::MatrixXd G(2, 2);
        G << 0.0, 0.0,  //
            0.0, 5.0;
        cases.push_back({"ZeroResidualGetsFullWeight", G,
                         Eigen::Vector2d(1.0, 0.0), 1e-12});
    }
    {
        // A single stored step leaves no freedom under sum(theta) = 1.
        Eigen::MatrixXd G(1, 1);
        G << 42.0;
        cases.push_back({"SingleHistoryEntryGetsWeightOne", G,
                         Eigen::VectorXd::Ones(1), 1e-12});
    }
    {
        // Identical residuals (f_0 == f_1) make the augmented system exactly
        // singular: every mixture has the same norm, so no minimizer is
        // distinguished. Expect the documented fallback, i.e. unit weight on
        // the newest step.
        Eigen::MatrixXd G(2, 2);
        G << 4.0, 4.0,  //
            4.0, 4.0;
        cases.push_back({"DuplicateResidualsFallBackToNewestStep", G,
                         Eigen::Vector2d(0.0, 1.0), 0.0});
    }
    {
        // f_1 is parallel to f_0 and differs from it only in the 9th
        // significant digit. The augmented system is singular to working
        // precision and the exact weights are O(1e9), i.e. pure cancellation
        // error. Expect the fallback rather than those weights.
        Eigen::MatrixXd F(2, 2);
        F << 1.0, 1.0 - 1e-9,  //
            0.0, 0.0;
        cases.push_back({"NearlyDuplicateResidualsFallBackToNewestStep",
                         gramOf(F), Eigen::Vector2d(0.0, 1.0), 0.0});
    }
    {
        // f_1 is parallel to f_0 and agrees with it to three digits. Here the
        // augmented system is still invertible and the solve returns the exact
        // weights [-999, 1000] - a lever arm that costs three digits of the
        // mixed iterate for no reliable gain, so the fallback must win.
        Eigen::MatrixXd F(2, 2);
        F << 1.0, 1.0 - 1e-3,  //
            0.0, 0.0;
        cases.push_back({"ExcessiveWeightsFallBackToNewestStep", gramOf(F),
                         Eigen::Vector2d(0.0, 1.0), 0.0});
    }
    {
        // All stored steps vanish, so there is nothing left to mix. The
        // normalization by the largest diagonal entry must not divide by zero.
        cases.push_back({"VanishingResidualsFallBackToNewestStep",
                         Eigen::MatrixXd::Zero(3, 3),
                         Eigen::Vector3d(0.0, 0.0, 1.0), 0.0});
    }
    {
        // The guards above must not fire on a healthy history: three decaying,
        // mutually orthogonal residuals, each a tenth of its predecessor. The
        // weights are the ones that zero out the model residual, concentrated
        // on the newest (smallest) step but with a genuine contribution from
        // the older ones - i.e. mixing actually takes place.
        Eigen::MatrixXd G(3, 3);
        G << 1.0, 0.0, 0.0,  //
            0.0, 1e-2, 0.0,  //
            0.0, 0.0, 1e-4;
        // theta_i proportional to 1/G_ii, normalized to sum one.
        Eigen::Vector3d const inverse_diagonal(1.0, 1e2, 1e4);
        cases.push_back({"DecayingIndependentResidualsAreMixed", G,
                         inverse_diagonal / inverse_diagonal.sum(), 1e-12});
    }

    return cases;
}

// The fixture stays in the unnamed namespace: its base
// testing::TestWithParam<AndersonWeightsCase> has internal linkage, so an
// externally visible fixture would trigger -Wsubobject-linkage.
struct NumLibAndersonAccelerationWeights
    : public ::testing::TestWithParam<AndersonWeightsCase>
{
};

TEST_P(NumLibAndersonAccelerationWeights, MatchExpectedWeights)
{
    auto const& c = GetParam();

    Eigen::VectorXd const theta =
        NumLib::detail::computeAndersonWeights(c.gram);

    ASSERT_EQ(theta.size(), c.expected_theta.size());
    for (Eigen::Index i = 0; i < theta.size(); ++i)
    {
        EXPECT_NEAR(theta(i), c.expected_theta(i), c.tolerance)
            << "component " << i;
    }
}

INSTANTIATE_TEST_SUITE_P(
    NumLib, NumLibAndersonAccelerationWeights,
    ::testing::ValuesIn(exactWeightsCases()),
    [](::testing::TestParamInfo<AndersonWeightsCase> const& info)
    { return info.param.name; });
}  // namespace

TEST(NumLibAndersonAcceleration, WeightsSumToOne)
{
    // Arbitrary symmetric positive definite Gram matrix (history size 3).
    Eigen::MatrixXd G(3, 3);
    G << 4.0, 1.0, 0.5,  //
        1.0, 3.0, 0.2,   //
        0.5, 0.2, 2.0;

    Eigen::VectorXd const theta = NumLib::detail::computeAndersonWeights(G);

    ASSERT_EQ(theta.size(), 3);
    EXPECT_NEAR(theta.sum(), 1.0, 1e-12);
}

TEST(NumLibAndersonAcceleration, MatchesBruteForceMinimizationOnRandomG)
{
    // Cross-check the closed-form solve against a brute-force scan over the
    // constrained simplex-like line theta_1 = 1 - theta_0 for a history of two,
    // i.e. minimize f(t) = ||t*f_0 + (1-t)*f_1||^2 = t^2*G00 + 2*t*(1-t)*G01 +
    // (1-t)^2*G11 over t directly and compare to the computed theta(0).
    Eigen::MatrixXd G(2, 2);
    G << 3.7, 1.3,  //
        1.3, 6.1;

    Eigen::VectorXd const theta = NumLib::detail::computeAndersonWeights(G);

    double const g00 = G(0, 0);
    double const g01 = G(0, 1);
    double const g11 = G(1, 1);
    // d/dt [t^2*g00 + 2t(1-t)*g01 + (1-t)^2*g11] = 0
    // => t*(g00 - 2*g01 + g11) = g11 - g01
    double const denom = g00 - 2 * g01 + g11;
    double const t_optimal = (g11 - g01) / denom;

    EXPECT_NEAR(theta(0), t_optimal, 1e-10);
    EXPECT_NEAR(theta(1), 1.0 - t_optimal, 1e-10);
}

TEST(NumLibAndersonAcceleration, ScaleInvarianceOfG)
{
    // theta must be unchanged when G is scaled by a positive constant, since
    // the norm-minimizing mixture direction is scale independent.
    Eigen::MatrixXd G(3, 3);
    G << 5.0, 2.0, 1.0,  //
        2.0, 4.0, 0.5,   //
        1.0, 0.5, 3.0;

    Eigen::VectorXd const theta_unscaled =
        NumLib::detail::computeAndersonWeights(G);
    Eigen::VectorXd const theta_scaled =
        NumLib::detail::computeAndersonWeights(1e6 * G);

    ASSERT_EQ(theta_unscaled.size(), theta_scaled.size());
    for (Eigen::Index i = 0; i < theta_unscaled.size(); ++i)
    {
        EXPECT_NEAR(theta_unscaled(i), theta_scaled(i), 1e-8);
    }
}

// The tests below drive the AndersonAcceleration class itself rather than the
// weight computation, i.e. the history buffer, its rotation and the
// incrementally maintained Gram matrix.
//
// Restricted to the serial build: constructing GlobalVector by length and
// reading entries back with get() is the serial EigenVector interface, whereas
// the PETSc vector distributes those indices over the ranks.
#ifndef USE_PETSC
namespace
{
//! A GlobalVector holding \c values, indexed from zero.
GlobalVector toGlobalVector(std::vector<double> const& values)
{
    GlobalVector v(static_cast<GlobalIndexType>(values.size()));
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        v.set(static_cast<GlobalIndexType>(i), values[i]);
    }
    return v;
}

//! Records one damped step x_old -> x_new and returns the (possibly mixed)
//! iterate the acceleration produces, mirroring how the Picard solver calls it.
std::vector<double> accelerated(NumLib::AndersonAcceleration& anderson,
                                std::vector<double> const& x_old,
                                std::vector<double> const& x_new)
{
    GlobalVector const x_old_vector = toGlobalVector(x_old);
    GlobalVector x_new_vector = toGlobalVector(x_new);

    anderson.accelerate(x_old_vector, x_new_vector);

    std::vector<double> result(x_new.size());
    for (std::size_t i = 0; i < result.size(); ++i)
    {
        result[i] = x_new_vector.get(static_cast<GlobalIndexType>(i));
    }
    return result;
}
}  // namespace

TEST(NumLibAndersonAccelerationHistory, FirstStepIsLeftUnchanged)
{
    NumLib::AndersonAcceleration anderson(2);

    // A single stored step carries unit weight by the sum-to-one constraint,
    // so the damped step the caller computed must come back untouched.
    auto const x_new = accelerated(anderson, {0.0, 0.0}, {1.0, 0.0});

    EXPECT_NEAR(x_new[0], 1.0, 1e-14);
    EXPECT_NEAR(x_new[1], 0.0, 1e-14);
}

TEST(NumLibAndersonAccelerationHistory, InertInstanceLeavesIterateUntouched)
{
    // A depth below min_mixing_depth stores nothing and mixes nothing.
    NumLib::AndersonAcceleration anderson(
        NumLib::AndersonAcceleration::min_mixing_depth - 1);

    accelerated(anderson, {0.0, 0.0}, {1.0, 0.0});
    auto const x_new = accelerated(anderson, {1.0, 0.0}, {1.0, 0.5});

    EXPECT_NEAR(x_new[0], 1.0, 1e-14);
    EXPECT_NEAR(x_new[1], 0.5, 1e-14);
}

TEST(NumLibAndersonAccelerationHistory, MixesTwoStepsWithExpectedWeights)
{
    NumLib::AndersonAcceleration anderson(2);

    // Steps f_0 = (1, 0) and f_1 = (0, 1/2) are orthogonal, so
    // G = diag(1, 1/4) and theta_i is proportional to 1/G_ii:
    // theta = (1, 4)/5 = (0.2, 0.8).
    accelerated(anderson, {0.0, 0.0}, {1.0, 0.0});
    auto const x_new = accelerated(anderson, {1.0, 0.0}, {1.0, 0.5});

    // 0.2*(x_0 + f_0) + 0.8*(x_1 + f_1) = 0.2*(1, 0) + 0.8*(1, 0.5).
    EXPECT_NEAR(x_new[0], 1.0, 1e-12);
    EXPECT_NEAR(x_new[1], 0.4, 1e-12);
}

TEST(NumLibAndersonAccelerationHistory, WindowRotatesAndGramFollowsIt)
{
    NumLib::AndersonAcceleration anderson(2);

    // Three steps into a window of two: the first is evicted, and the mixture
    // must be the one of steps two and three alone. A Gram matrix that was not
    // shifted along with the buffer would mix the stale first step's norm in
    // and produce different weights.
    accelerated(anderson, {0.0, 0.0}, {1.0, 0.0});
    accelerated(anderson, {10.0, 0.0}, {10.0, 0.5});
    auto const x_new = accelerated(anderson, {20.0, 0.0}, {20.125, 0.0});

    // f_1 = (0, 1/2), f_2 = (1/8, 0) are orthogonal, so G = diag(1/4, 1/64)
    // and theta is proportional to (4, 64), i.e. (4, 64)/68.
    double const theta_1 = 4.0 / 68.0;
    double const theta_2 = 64.0 / 68.0;
    EXPECT_NEAR(x_new[0], theta_1 * 10.0 + theta_2 * 20.125, 1e-12);
    EXPECT_NEAR(x_new[1], theta_1 * 0.5, 1e-12);
}

TEST(NumLibAndersonAccelerationHistory, DropLastStepUndoesTheRecordedStep)
{
    NumLib::AndersonAcceleration with_dropped_step(2);
    NumLib::AndersonAcceleration reference(2);

    accelerated(with_dropped_step, {0.0, 0.0}, {1.0, 0.0});
    accelerated(reference, {0.0, 0.0}, {1.0, 0.0});

    // A repeated iteration: the step is recorded, then dropped again.
    accelerated(with_dropped_step, {5.0, 5.0}, {6.0, 4.0});
    with_dropped_step.dropLastStep();

    auto const dropped = accelerated(with_dropped_step, {1.0, 0.0}, {1.0, 0.5});
    auto const expected = accelerated(reference, {1.0, 0.0}, {1.0, 0.5});

    ASSERT_EQ(dropped.size(), expected.size());
    for (std::size_t i = 0; i < dropped.size(); ++i)
    {
        EXPECT_NEAR(dropped[i], expected[i], 1e-14);
    }
}
#endif  // USE_PETSC
