/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <Eigen/Core>
#include <random>

#include "NumLib/NumericalDifferentiation.h"
#include "Tests/TestTools.h"

namespace
{
// only scalar arguments
Eigen::Vector2d fNumDiffScalar(double const t, double const x, double const y,
                               double const z)
{
    return {2 * std::sin(3 * x - 5 * y + 7 * z - 11 * t),
            3 * std::cos(5 * x - 7 * y + 11 * z - 2 * t)};
}

Eigen::Vector2d dfNumDiffScalar_dt(double const t, double const x,
                                   double const y, double const z)
{
    return {-22 * std::cos(3 * x - 5 * y + 7 * z - 11 * t),
            6 * std::sin(5 * x - 7 * y + 11 * z - 2 * t)};
}

Eigen::Vector2d dfNumDiffScalar_dx(double const t, double const x,
                                   double const y, double const z)
{
    return {6 * std::cos(3 * x - 5 * y + 7 * z - 11 * t),
            -15 * std::sin(5 * x - 7 * y + 11 * z - 2 * t)};
}

Eigen::Vector2d dfNumDiffScalar_dy(double const t, double const x,
                                   double const y, double const z)
{
    return {-10 * std::cos(3 * x - 5 * y + 7 * z - 11 * t),
            21 * std::sin(5 * x - 7 * y + 11 * z - 2 * t)};
}

Eigen::Vector2d dfNumDiffScalar_dz(double const t, double const x,
                                   double const y, double const z)
{
    return {14 * std::cos(3 * x - 5 * y + 7 * z - 11 * t),
            -33 * std::sin(5 * x - 7 * y + 11 * z - 2 * t)};
}

// also vectorial arguments
Eigen::Vector2d fNumDiffVectorial(double const t, Eigen::Vector2d const& xy,
                                  double const z)
{
    return fNumDiffScalar(t, xy[0], xy[1], z);
}

struct NumDiffXY
{
    NumDiffXY(double x, double y) : x{x}, y{y} {}
    explicit NumDiffXY(Eigen::Vector2d const& xy) : x{xy[0]}, y{xy[1]} {}

    NumDiffXY operator-(NumDiffXY const& other) const
    {
        return {x - other.x, y - other.y};
    }

    NumDiffXY operator/(double a) const { return {x / a, y / a}; }

    double x;
    double y;
};

// returns a custom struct
NumDiffXY fNumDiffCustom(double const t, Eigen::Vector3d const& xyz)
{
    return NumDiffXY{fNumDiffScalar(t, xyz[0], xyz[1], xyz[2])};
}
}  // namespace

// Test fixture.
template <typename T>
class NumLibNumericalDifferentiationTest : public ::testing::Test
{
protected:
    static typename T::Derivative createDerivative()
    {
        return {NumLib::RelativeEpsilon{rel_eps},
                NumLib::MinimumPerturbation{min_pert}};
    }

    double randomValue()
    {
        auto const sign = sign_distribution(engine) == 0 ? -1 : 1;
        return sign * magnitude_distribution(engine);
    }

    void runStandardVectorialCase(typename T::Derivative const& deriv,
                                  double const t,
                                  double const x,
                                  double const y,
                                  double const z)
    {
        auto const f = fNumDiffVectorial(t, Eigen::Vector2d{x, y}, z);
        auto const [dfdt, dfdxy, dfdz] =
            deriv(&fNumDiffVectorial, t, Eigen::Vector2d{x, y}, z);

        static_assert(std::is_same_v<decltype(f), decltype(dfdt)>);
        using VectorialReturn = std::array<std::remove_const_t<decltype(f)>, 2>;
        static_assert(std::is_same_v<const VectorialReturn, decltype(dfdxy)>);
        static_assert(std::is_same_v<decltype(f), decltype(dfdz)>);

        constexpr auto abstol = T::abstol;
        EXPECT_PRED_FORMAT3(Tests::EigenIsNear{}, dfdt,
                            dfNumDiffScalar_dt(t, x, y, z), abstol);
        EXPECT_PRED_FORMAT3(Tests::EigenIsNear{}, dfdxy[0],
                            dfNumDiffScalar_dx(t, x, y, z), abstol);
        EXPECT_PRED_FORMAT3(Tests::EigenIsNear{}, dfdxy[1],
                            dfNumDiffScalar_dy(t, x, y, z), abstol);
        EXPECT_PRED_FORMAT3(Tests::EigenIsNear{}, dfdz,
                            dfNumDiffScalar_dz(t, x, y, z), abstol);
    }

    static constexpr double rel_eps = 1e-8;
    static constexpr double min_pert = 1e-8;
    static constexpr int num_repetitions = 100;

private:
    std::random_device device;
    std::default_random_engine engine{device()};

    // range chosen such that the produced values lie in the "normal"
    // perturbation range of the tested differentiation algorithms
    std::uniform_real_distribution<double> magnitude_distribution{
        min_pert / rel_eps, 10 * min_pert / rel_eps};

    std::uniform_int_distribution<int> sign_distribution{0, 1};
};

struct NumLibNumericalDifferentiationTestCD
{
    using Derivative =
        NumLib::NumericalDerivative<NumLib::CentralDifferencesStrategy>;

    static constexpr double abstol = 9e-6;
};

struct NumLibNumericalDifferentiationTestFD
{
    using Derivative =
        NumLib::NumericalDerivative<NumLib::ForwardDifferencesStrategy>;

    static constexpr double abstol = 2e-5;
};

using NumLibNumericalDifferentiationTestTypes =
    ::testing::Types<NumLibNumericalDifferentiationTestCD,
                     NumLibNumericalDifferentiationTestFD>;

TYPED_TEST_SUITE(NumLibNumericalDifferentiationTest,
                 NumLibNumericalDifferentiationTestTypes);

TYPED_TEST(NumLibNumericalDifferentiationTest, AllScalarWithLambda)
{
    auto const deriv = this->createDerivative();

    for (int repetition = 0; repetition < this->num_repetitions; ++repetition)
    {
        double const t = this->randomValue();
        double const x = this->randomValue();
        double const y = this->randomValue();
        double const z = this->randomValue();

        auto const f = fNumDiffScalar(t, x, y, z);
        auto const [dfdt, dfdx, dfdy, dfdz] = deriv(
            [](auto... args) { return fNumDiffScalar(args...); }, t, x, y, z);

        static_assert(std::is_same_v<decltype(f), decltype(dfdt)>);
        static_assert(std::is_same_v<decltype(f), decltype(dfdx)>);
        static_assert(std::is_same_v<decltype(f), decltype(dfdy)>);
        static_assert(std::is_same_v<decltype(f), decltype(dfdz)>);

        constexpr auto abstol = TypeParam::abstol;
        EXPECT_PRED_FORMAT3(Tests::EigenIsNear{}, dfdt,
                            dfNumDiffScalar_dt(t, x, y, z), abstol);
        EXPECT_PRED_FORMAT3(Tests::EigenIsNear{}, dfdx,
                            dfNumDiffScalar_dx(t, x, y, z), abstol);
        EXPECT_PRED_FORMAT3(Tests::EigenIsNear{}, dfdy,
                            dfNumDiffScalar_dy(t, x, y, z), abstol);
        EXPECT_PRED_FORMAT3(Tests::EigenIsNear{}, dfdz,
                            dfNumDiffScalar_dz(t, x, y, z), abstol);
    }
}

TYPED_TEST(NumLibNumericalDifferentiationTest, Vectorial)
{
    auto const deriv = this->createDerivative();

    for (int repetition = 0; repetition < this->num_repetitions; ++repetition)
    {
        double const t = this->randomValue();
        double const x = this->randomValue();
        double const y = this->randomValue();
        double const z = this->randomValue();

        this->runStandardVectorialCase(deriv, t, x, y, z);
    }
}

TYPED_TEST(NumLibNumericalDifferentiationTest, CustomStruct)
{
    auto const deriv = this->createDerivative();

    for (int repetition = 0; repetition < this->num_repetitions; ++repetition)
    {
        double const t = this->randomValue();
        double const x = this->randomValue();
        double const y = this->randomValue();
        double const z = this->randomValue();

        auto const f = fNumDiffCustom(t, Eigen::Vector3d{x, y, z});
        auto const [dfdt, grad_f] =
            deriv(&fNumDiffCustom, t, Eigen::Vector3d{x, y, z});

        static_assert(std::is_same_v<decltype(f), decltype(dfdt)>);
        using VectorialReturn = std::array<std::remove_const_t<decltype(f)>, 3>;
        static_assert(std::is_same_v<const VectorialReturn, decltype(grad_f)>);

        constexpr auto abstol = TypeParam::abstol;
        {
            auto const expected = dfNumDiffScalar_dt(t, x, y, z);
            EXPECT_NEAR(expected[0], dfdt.x, abstol);
            EXPECT_NEAR(expected[1], dfdt.y, abstol);
        }

        {
            auto const expected = dfNumDiffScalar_dx(t, x, y, z);
            EXPECT_NEAR(expected[0], grad_f[0].x, abstol);
            EXPECT_NEAR(expected[1], grad_f[0].y, abstol);
        }

        {
            auto const expected = dfNumDiffScalar_dy(t, x, y, z);
            EXPECT_NEAR(expected[0], grad_f[1].x, abstol);
            EXPECT_NEAR(expected[1], grad_f[1].y, abstol);
        }

        {
            auto const expected = dfNumDiffScalar_dz(t, x, y, z);
            EXPECT_NEAR(expected[0], grad_f[2].x, abstol);
            EXPECT_NEAR(expected[1], grad_f[2].y, abstol);
        }
    }
}

TYPED_TEST(NumLibNumericalDifferentiationTest, NegativePerturbation)
{
    typename TypeParam::Derivative deriv{
        NumLib::RelativeEpsilon{-this->rel_eps},
        NumLib::MinimumPerturbation{-this->min_pert}};

    for (int repetition = 0; repetition < this->num_repetitions; ++repetition)
    {
        double const t = this->randomValue();
        double const x = this->randomValue();
        double const y = this->randomValue();
        double const z = this->randomValue();

        this->runStandardVectorialCase(deriv, t, x, y, z);
    }
}

TYPED_TEST(NumLibNumericalDifferentiationTest, CloseToZero)
{
    auto const run = [this](double const rel_eps, double const min_pert)
    {
        typename TypeParam::Derivative deriv{
            NumLib::RelativeEpsilon{rel_eps},
            NumLib::MinimumPerturbation{min_pert}};

        double t = this->randomValue();
        double x = this->randomValue();
        double y = this->randomValue();
        double z = this->randomValue();

        std::array<std::reference_wrapper<double>, 4> all_inputs{t, x, y, z};

        for (auto const& input : all_inputs)
        {
            double const old_value = input;

            double const a = min_pert / rel_eps;

            for (double const new_value : {0.0, 0.5 * min_pert, min_pert,
                                           2 * min_pert, 0.5 * a, a, 2 * a})
            {
                input.get() = new_value;

                this->runStandardVectorialCase(deriv, t, x, y, z);
            }

            // reset input
            input.get() = old_value;
        }
    };

    for (int repetition = 0; repetition < this->num_repetitions; ++repetition)
    {
        for (double const rel_eps : {this->rel_eps, -this->rel_eps})
        {
            for (double const min_pert : {this->min_pert, -this->min_pert})
            {
                run(rel_eps, min_pert);
            }
        }
    }
}
