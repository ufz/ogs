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

#include "MaterialLib/MPL/Properties/Curve.h"

TEST(MaterialPropertyLib, Curve)
{
    MathLib::PiecewiseLinearInterpolation curves{
        {300, 307},  // coords
        {300, 500}   // values
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
        {0, 100},   // coords
        {300, 500}  // values
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
        {0, 100},   // coords
        {300, 500}  // values
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
