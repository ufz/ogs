/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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
    variable_array.temperature = 303.15;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    ASSERT_NEAR(
        std::get<double>(linear_property.value(variable_array, pos, time, dt)),
        y_ref * (1 + m * (variable_array.temperature - x_ref)),
        1.e-10);
    ASSERT_EQ(std::get<double>(linear_property.dValue(
                  variable_array,
                  MaterialPropertyLib::Variable::liquid_phase_pressure, pos,
                  time, dt)),
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

TEST(MaterialPropertyLib, LinearWithTime)
{
    double const y_ref = 1.0;
    double const m = 1.0;
    double const x_ref = 0.0;

    MaterialPropertyLib::IndependentVariable const iv{"t", x_ref, m};

    std::vector<MaterialPropertyLib::IndependentVariable> ivs{iv};
    MaterialPropertyLib::Linear linear_property{"linear", y_ref, ivs};

    MaterialPropertyLib::VariableArray variable_array;

    ParameterLib::SpatialPosition const pos;
    double const time = 10;
    double const dt = std::numeric_limits<double>::quiet_NaN();
    ASSERT_NEAR(
        std::get<double>(linear_property.value(variable_array, pos, time, dt)),
        y_ref * (1 + m * (time - x_ref)),
        1.e-10);
    ASSERT_EQ(std::get<double>(linear_property.dValue(
                  variable_array,
                  MaterialPropertyLib::Variable::liquid_phase_pressure, pos,
                  time, dt)),
              0.0);
}

TEST(MaterialPropertyLib, LinearWithPos)
{
    double const y_ref = 1.0;
    double const m = 1.0;
    double const x_ref = 0.0;

    MaterialPropertyLib::IndependentVariable const iv_x{"x", x_ref, m};
    MaterialPropertyLib::IndependentVariable const iv_y{"y", x_ref, m};
    MaterialPropertyLib::IndependentVariable const iv_z{"z", x_ref, m};

    std::vector<MaterialPropertyLib::IndependentVariable> ivs_x{iv_x};
    MaterialPropertyLib::Linear linear_property_x{"linear", y_ref, ivs_x};
    std::vector<MaterialPropertyLib::IndependentVariable> ivs_y{iv_y};
    MaterialPropertyLib::Linear linear_property_y{"linear", y_ref, ivs_y};
    std::vector<MaterialPropertyLib::IndependentVariable> ivs_z{iv_z};
    MaterialPropertyLib::Linear linear_property_z{"linear", y_ref, ivs_z};

    MaterialPropertyLib::VariableArray variable_array;

    std::array<double, 3> coords = {50, 20, -10};
    ParameterLib::SpatialPosition const pos =
        ParameterLib::SpatialPosition{{0}, {0}, MathLib::Point3d{coords}};
    double const time = 10;
    double const dt = std::numeric_limits<double>::quiet_NaN();

    ASSERT_NEAR(std::get<double>(
                    linear_property_x.value(variable_array, pos, time, dt)),
                y_ref * (1 + m * (coords[0] - x_ref)),
                1.e-10);
    ASSERT_EQ(std::get<double>(linear_property_x.dValue(
                  variable_array,
                  MaterialPropertyLib::Variable::liquid_phase_pressure, pos,
                  time, dt)),
              0.0);
    ASSERT_NEAR(std::get<double>(
                    linear_property_y.value(variable_array, pos, time, dt)),
                y_ref * (1 + m * (coords[1] - x_ref)),
                1.e-10);
    ASSERT_EQ(std::get<double>(linear_property_y.dValue(
                  variable_array,
                  MaterialPropertyLib::Variable::liquid_phase_pressure, pos,
                  time, dt)),
              0.0);
    ASSERT_NEAR(std::get<double>(
                    linear_property_z.value(variable_array, pos, time, dt)),
                y_ref * (1 + m * (coords[2] - x_ref)),
                1.e-10);
    ASSERT_EQ(std::get<double>(linear_property_z.dValue(
                  variable_array,
                  MaterialPropertyLib::Variable::liquid_phase_pressure, pos,
                  time, dt)),
              0.0);
}