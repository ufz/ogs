/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <limits>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/MPL/Properties/Function.h"
#include "MathLib/KelvinVector.h"
#include "MathLib/VectorizedTensor.h"

namespace MPL = MaterialPropertyLib;

struct MPLFunction : public ::testing::Test
{
    static double constexpr nan = std::numeric_limits<double>::quiet_NaN();
    MPL::VariableArray vars;
};

TEST_F(MPLFunction, ScalarScalar)
{
    vars.temperature = 2.;

    MPL::Property const& f = MPL::Function{
        "test_function", {"temperature"}, {{"temperature", {"1"}}}};
    ASSERT_EQ(2., f.value<double>(vars, {}, nan, nan));
    ASSERT_EQ(1.,
              f.dValue<double>(vars, MPL::Variable::temperature, {}, nan, nan));
}

TEST_F(MPLFunction, ScalarVector)
{
    vars.temperature = 2.;

    MPL::Property const& f =
        MPL::Function{"test_function",
                      {"temperature", "temperature^2"},
                      {{"temperature", {"1", "2*temperature"}}}};
    ASSERT_EQ((Eigen::Vector2d{2., 4.}),
              (f.value<Eigen::Vector2d>(vars, {}, nan, nan)));
    ASSERT_EQ((Eigen::Vector2d{1., 4.}),
              f.dValue<Eigen::Vector2d>(
                  vars, MPL::Variable::temperature, {}, nan, nan));
}

TEST_F(MPLFunction, ScalarUninitialized)
{
    // The vars.temperature = is not initialized.

    MPL::Property const& f = MPL::Function{
        "test_function", {"temperature"}, {{"temperature", {"1"}}}};
    ASSERT_ANY_THROW(f.value<double>(vars, {}, nan, nan));
    ASSERT_ANY_THROW(
        f.dValue<double>(vars, MPL::Variable::temperature, {}, nan, nan));
}
