/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include "MaterialLib/MPL/Properties/Linear.h"

TEST(MaterialPropertyLib, Linear)
{
    double const y_ref = 1.0;
    double const m = 1.0;
    double const x_ref = 293.15;
    MaterialPropertyLib::IndependentVariable const iv{
        MaterialPropertyLib::Variable::temperature, x_ref, m};

    std::vector<MaterialPropertyLib::IndependentVariable> ivs{iv};
    MaterialPropertyLib::Linear linear_property{"linear", y_ref, ivs};

    MaterialPropertyLib::VariableArray variable_array;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::temperature)] = 303.15;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    ASSERT_NEAR(
        std::get<double>(linear_property.value(variable_array, pos, time, dt)),
        y_ref * (1 + m * (std::get<double>(variable_array[static_cast<int>(
                              MaterialPropertyLib::Variable::temperature)]) -
                          x_ref)),
        1.e-10);
    ASSERT_EQ(std::get<double>(linear_property.dValue(
                  variable_array, MaterialPropertyLib::Variable::phase_pressure,
                  pos, time, dt)),
              0.0);
    ASSERT_NEAR(std::get<double>(linear_property.dValue(
                    variable_array, MaterialPropertyLib::Variable::temperature,
                    pos, time, dt)),
                y_ref * m, 1.e-16);
    ASSERT_EQ(std::get<double>(linear_property.d2Value(
                  variable_array, MaterialPropertyLib::Variable::temperature,
                  MaterialPropertyLib::Variable::temperature, pos, time, dt)),
              0.0);
}
