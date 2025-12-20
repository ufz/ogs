// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include "MaterialLib/MPL/Properties/Curve.h"

TEST(MaterialPropertyLib, Curve)
{
    MathLib::PiecewiseLinearInterpolation curves{
        std::vector<int>({300, 307}),    // coords
        std::vector<double>({300, 500})  // values
    };

    MaterialPropertyLib::VariableArray variable_array;
    variable_array.temperature = 303.5;
    MaterialPropertyLib::Property const& curve_property =
        MaterialPropertyLib::Curve{
            "test_curve", MaterialPropertyLib::Variable::temperature, curves};
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    ASSERT_EQ(
        std::get<double>(curve_property.value(variable_array, pos, time, dt)),
        400);
    ASSERT_EQ(std::get<double>(curve_property.dValue(
                  variable_array,
                  MaterialPropertyLib::Variable::temperature,
                  pos,
                  time,
                  dt)),
              200. / 7.);
}

TEST(MaterialPropertyLib, CurveTime)
{
    MathLib::PiecewiseLinearInterpolation curves{
        std::vector<double>({0, 100}),   // coords
        std::vector<double>({300, 500})  // values
    };

    MaterialPropertyLib::VariableArray variable_array;
    MaterialPropertyLib::Property const& curve_property =
        MaterialPropertyLib::Curve{"test_curve", "t", curves};
    ParameterLib::SpatialPosition const pos;
    double const time = 50;
    double const dt = std::numeric_limits<double>::quiet_NaN();
    ASSERT_EQ(
        std::get<double>(curve_property.value(variable_array, pos, time, dt)),
        400);
    ASSERT_EQ(std::get<double>(curve_property.dValue(
                  variable_array,
                  MaterialPropertyLib::Variable::temperature,
                  pos,
                  time,
                  dt)),
              0);
}

TEST(MaterialPropertyLib, CurvePos)
{
    MathLib::PiecewiseLinearInterpolation curves{
        std::vector<int>({0, 100}),      // coords
        std::vector<double>({300, 500})  // values
    };

    MaterialPropertyLib::VariableArray variable_array;
    MaterialPropertyLib::Property const& curve_property_x =
        MaterialPropertyLib::Curve{"test_curve", "x", curves};
    MaterialPropertyLib::Property const& curve_property_y =
        MaterialPropertyLib::Curve{"test_curve", "y", curves};
    MaterialPropertyLib::Property const& curve_property_z =
        MaterialPropertyLib::Curve{"test_curve", "z", curves};
    std::array<double, 3> coords = {50, 25, 75};
    ParameterLib::SpatialPosition const pos =
        ParameterLib::SpatialPosition{{0}, {0}, MathLib::Point3d{coords}};
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    ASSERT_EQ(
        std::get<double>(curve_property_x.value(variable_array, pos, time, dt)),
        400);
    ASSERT_EQ(std::get<double>(curve_property_x.dValue(
                  variable_array,
                  MaterialPropertyLib::Variable::temperature,
                  pos,
                  time,
                  dt)),
              0);
    ASSERT_EQ(
        std::get<double>(curve_property_y.value(variable_array, pos, time, dt)),
        350);
    ASSERT_EQ(std::get<double>(curve_property_y.dValue(
                  variable_array,
                  MaterialPropertyLib::Variable::temperature,
                  pos,
                  time,
                  dt)),
              0);
    ASSERT_EQ(
        std::get<double>(curve_property_z.value(variable_array, pos, time, dt)),
        450);
    ASSERT_EQ(std::get<double>(curve_property_z.dValue(
                  variable_array,
                  MaterialPropertyLib::Variable::temperature,
                  pos,
                  time,
                  dt)),
              0);
}
