/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <gtest/gtest.h>

#include "MaterialLib/MPL/Utils/SigmoidFunction.h"

namespace MaterialPropertyLib
{
class SigmoidFunctionTest : public ::testing::Test
{
protected:
    SigmoidFunction sigmoid_function_with_residual_saturation{1.0, 0.0, 0.1};
    SigmoidFunction sigmoid_function_without_residual_saturation{1.0, 0.0, 0.0};
};

TEST_F(SigmoidFunctionTest, ValueTest)
{
    EXPECT_NEAR(sigmoid_function_with_residual_saturation.value(0.0),
                (1.0 - 0.1) / (1.0 + std::exp(0.0)), 1e-6);
    EXPECT_NEAR(sigmoid_function_with_residual_saturation.value(-1.0),
                (1.0 - 0.1) / (1.0 + std::exp(-1.0)), 1e-6);
    EXPECT_NEAR(sigmoid_function_with_residual_saturation.value(1.0),
                (1.0 - 0.1) / (1.0 + std::exp(1.0)), 1e-6);

    EXPECT_NEAR(sigmoid_function_without_residual_saturation.value(0.0),
                1.0 / (1.0 + std::exp(0.0)), 1e-6);
    EXPECT_NEAR(sigmoid_function_without_residual_saturation.value(-1.0),
                1.0 / (1.0 + std::exp(-1.0)), 1e-6);
    EXPECT_NEAR(sigmoid_function_without_residual_saturation.value(1.0),
                1.0 / (1.0 + std::exp(1.0)), 1e-6);
}

TEST_F(SigmoidFunctionTest, DValueTest)
{
    double scale_factor = (1.0 - 0.1);
    EXPECT_NEAR(
        sigmoid_function_with_residual_saturation.dValue(0.0),
        scale_factor * (-std::exp(0.0) / std::pow(1.0 + std::exp(0.0), 2)),
        1e-6);
    EXPECT_NEAR(
        sigmoid_function_with_residual_saturation.dValue(1.0),
        scale_factor * (-std::exp(1.0) / std::pow(1.0 + std::exp(1.0), 2)),
        1e-6);
    EXPECT_NEAR(
        sigmoid_function_with_residual_saturation.dValue(-1.0),
        scale_factor * (-std::exp(-1.0) / std::pow(1.0 + std::exp(-1.0), 2)),
        1e-6);

    EXPECT_NEAR(sigmoid_function_without_residual_saturation.dValue(0.0),
                -std::exp(0.0) / std::pow(1.0 + std::exp(0.0), 2), 1e-6);
    EXPECT_NEAR(sigmoid_function_without_residual_saturation.dValue(1.0),
                -std::exp(1.0) / std::pow(1.0 + std::exp(1.0), 2), 1e-6);
    EXPECT_NEAR(sigmoid_function_without_residual_saturation.dValue(-1.0),
                -std::exp(-1.0) / std::pow(1.0 + std::exp(-1.0), 2), 1e-6);
}

TEST_F(SigmoidFunctionTest, D2ValueTest)
{
    double const k = 1.0;
    double const scale_factor =
        1.0 - 0.1;  // for sigmoid_function_with_residual_saturation

    EXPECT_NEAR(sigmoid_function_with_residual_saturation.d2Value(0.0),
                scale_factor * k * k * std::exp(0.0) * (std::exp(0.0) - 1.0) /
                    std::pow(1.0 + std::exp(0.0), 3),
                1e-6);
    EXPECT_NEAR(sigmoid_function_with_residual_saturation.d2Value(1.0),
                scale_factor * k * k * std::exp(1.0) * (std::exp(1.0) - 1.0) /
                    std::pow(1.0 + std::exp(1.0), 3),
                1e-6);
    EXPECT_NEAR(sigmoid_function_with_residual_saturation.d2Value(-1.0),
                scale_factor * k * k * std::exp(-1.0) * (std::exp(-1.0) - 1.0) /
                    std::pow(1.0 + std::exp(-1.0), 3),
                1e-6);

    EXPECT_NEAR(sigmoid_function_without_residual_saturation.d2Value(0.0),
                k * k * std::exp(0.0) * (std::exp(0.0) - 1.0) /
                    std::pow(1.0 + std::exp(0.0), 3),
                1e-6);
    EXPECT_NEAR(sigmoid_function_without_residual_saturation.d2Value(1.0),
                k * k * std::exp(1.0) * (std::exp(1.0) - 1.0) /
                    std::pow(1.0 + std::exp(1.0), 3),
                1e-6);
    EXPECT_NEAR(sigmoid_function_without_residual_saturation.d2Value(-1.0),
                k * k * std::exp(-1.0) * (std::exp(-1.0) - 1.0) /
                    std::pow(1.0 + std::exp(-1.0), 3),
                1e-6);
}

}  // namespace MaterialPropertyLib
