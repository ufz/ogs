// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <cmath>
#include <limits>

#include "MaterialLib/MPL/Properties/Sigmoid.h"
#include "TestMPL.h"
#include "Tests/TestTools.h"

// Helpers shared by all tests. The properties used here keep steepness k = 1
// and midpoint X_c = 0, so that with lower_bound = 0 and
// scale = upper_bound the closed-form expressions (monotonically increasing)
// are
//   value   = scale / (1 + exp(-X))
//   dValue  = scale * exp(-X) / (1 + exp(-X))^2
//   d2Value = scale * exp(-X) * (exp(-X) - 1) / (1 + exp(-X))^3
namespace
{
constexpr double k = 1.0;
constexpr double X_c = 0.0;
constexpr double S_r = 0.1;

MaterialPropertyLib::Sigmoid makeSigmoid(double const lower_bound,
                                         double const upper_bound)
{
    return MaterialPropertyLib::Sigmoid("frozen_liquid_saturation", k, X_c,
                                        lower_bound, upper_bound,
                                        MPL::Variable::temperature);
}

double valueAt(MPL::Property const& s, double const X)
{
    ParameterLib::SpatialPosition const pos;
    double const nan = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;
    vars.temperature = X;
    return s.value<double>(vars, pos, nan, nan);
}

double dValueAt(MPL::Property const& s, double const X,
                MPL::Variable const variable)
{
    ParameterLib::SpatialPosition const pos;
    double const nan = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;
    vars.temperature = X;
    return s.dValue<double>(vars, variable, pos, nan, nan);
}

double d2ValueAt(MPL::Property const& s, double const X)
{
    ParameterLib::SpatialPosition const pos;
    double const nan = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;
    vars.temperature = X;
    return s.d2Value<double>(vars, MPL::Variable::temperature,
                             MPL::Variable::temperature, pos, nan, nan);
}
}  // namespace

// Checks against the closed-form formulae for value, dValue and d2Value,
// sampled at several points. Two properties are used: one with a residual
// saturation (upper_bound = 1 - S_r) and one without (upper_bound = 1). The
// sample set includes X = X_c = 0, so these also cover the midpoint
// (value = 0.5, dValue = +k/4, d2Value = 0 for the property without residual
// saturation).
TEST(MaterialPropertyLib, Sigmoid_value_closedForm)
{
    auto const with_residual = makeSigmoid(0.0, 1.0 - S_r);
    auto const without_residual = makeSigmoid(0.0, 1.0);

    for (double const X : {0.0, 1.0, -1.0})
    {
        EXPECT_NEAR(valueAt(with_residual, X),
                    (1.0 - S_r) / (1.0 + std::exp(-X)), 1e-6);
        EXPECT_NEAR(valueAt(without_residual, X), 1.0 / (1.0 + std::exp(-X)),
                    1e-6);
    }
}

TEST(MaterialPropertyLib, Sigmoid_dValue_closedForm)
{
    auto const with_residual = makeSigmoid(0.0, 1.0 - S_r);
    auto const without_residual = makeSigmoid(0.0, 1.0);

    for (double const X : {0.0, 1.0, -1.0})
    {
        EXPECT_NEAR(
            dValueAt(with_residual, X, MPL::Variable::temperature),
            (1.0 - S_r) * (std::exp(-X) / std::pow(1.0 + std::exp(-X), 2)),
            1e-6);
        EXPECT_NEAR(dValueAt(without_residual, X, MPL::Variable::temperature),
                    std::exp(-X) / std::pow(1.0 + std::exp(-X), 2), 1e-6);
    }
}

TEST(MaterialPropertyLib, Sigmoid_d2Value_closedForm)
{
    auto const with_residual = makeSigmoid(0.0, 1.0 - S_r);
    auto const without_residual = makeSigmoid(0.0, 1.0);

    for (double const X : {0.0, 1.0, -1.0})
    {
        EXPECT_NEAR(d2ValueAt(with_residual, X),
                    (1.0 - S_r) * std::exp(-X) * (std::exp(-X) - 1.0) /
                        std::pow(1.0 + std::exp(-X), 3),
                    1e-6);
        EXPECT_NEAR(d2ValueAt(without_residual, X),
                    std::exp(-X) * (std::exp(-X) - 1.0) /
                        std::pow(1.0 + std::exp(-X), 3),
                    1e-6);
    }
}

TEST(MaterialPropertyLib, Sigmoid_value_extremes)
{
    auto const sigmoid = makeSigmoid(0.0, 1.0);

    // Far below characteristic value (X << X_c): 1 / (1 + exp(100)) ~= 0
    EXPECT_LT(valueAt(sigmoid, X_c - 100.0), 0.01);
    // Far above characteristic value (X >> X_c): 1 / (1 + exp(-100)) ~= 1
    EXPECT_GT(valueAt(sigmoid, X_c + 100.0), 0.99);
}

TEST(MaterialPropertyLib, Sigmoid_value_overflow_guard)
{
    // With k = 1 and X_c = 0 the internal exponent is -(X - X_c) = -X. Choosing
    // X far below the midpoint drives the exponent past the std::exp overflow
    // threshold (max_exp_argument ~= 709.78), exercising the guard in value().
    // The result must stay finite and collapse exactly onto the lower bound.
    auto const sigmoid = makeSigmoid(0.2, 0.8);

    double const value = valueAt(sigmoid, X_c - 1000.0);
    EXPECT_TRUE(std::isfinite(value));
    EXPECT_EQ(value, 0.2);
}

TEST(MaterialPropertyLib, Sigmoid_with_custom_bounds)
{
    // At the midpoint with custom bounds:
    // (0.8 - 0.2) / (1 + 1) + 0.2 = 0.5
    auto const sigmoid = makeSigmoid(0.2, 0.8);

    EXPECT_NEAR(valueAt(sigmoid, X_c), 0.5, 1e-10);

    // At the midpoint the derivative equals steepness * range / 4, with
    // range = upper_bound - lower_bound = 0.6.
    EXPECT_NEAR(dValueAt(sigmoid, X_c, MPL::Variable::temperature),
                k * (0.8 - 0.2) / 4.0, 1e-10);
}

TEST(MaterialPropertyLib, Sigmoid_derivative_zero_for_unrelated_variable)
{
    auto const sigmoid = makeSigmoid(0.0, 1.0);

    // Derivative with respect to capillary_pressure (not the independent
    // variable) must be zero.
    EXPECT_EQ(dValueAt(sigmoid, X_c, MPL::Variable::capillary_pressure), 0.0);
}
