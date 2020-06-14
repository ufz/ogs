/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include <cmath>

#include "MaterialLib/MPL/Properties/Exponential.h"

namespace MPL = MaterialPropertyLib;

TEST(MaterialPropertyLib, Exponential)
{
    double const y_ref = 1e-3;
    double const reference_condition = 20.0;
    double const factor = 1 / 75.0;
    MPL::ExponentData const exp_data{MPL::Variable::temperature,
                                     reference_condition, factor};
    MPL::Property const& p = MPL::Exponential{"exponential", y_ref, exp_data};

    double const T = 20.;
    MPL::VariableArray variable_array;
    variable_array[static_cast<int>(MPL::Variable::temperature)] = T;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    ASSERT_NEAR(p.template value<double>(variable_array, pos, time, dt),
                y_ref * (std::exp(factor * (T - reference_condition))),
                1.e-10);
    ASSERT_EQ(p.template dValue<double>(
                  variable_array, MPL::Variable::phase_pressure, pos, time, dt),
              0.0);
    ASSERT_NEAR(p.template dValue<double>(
                    variable_array, MPL::Variable::temperature, pos, time, dt),
                y_ref * factor * std::exp(factor * (T - reference_condition)),
                1.e-16);
    ASSERT_NEAR(
        p.template d2Value<double>(variable_array, MPL::Variable::temperature,
                                   MPL::Variable::temperature, pos, time, dt),
        y_ref * std::pow(factor, 2) *
            std::exp(factor * (T - reference_condition)),
        1.e-16);
}
