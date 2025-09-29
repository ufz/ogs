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

#include <range/v3/view/enumerate.hpp>

#include "MaterialLib/MPL/Properties/Linear.h"

struct MaterialPropertyLibLinearTestParams
{
    std::optional<double> min;
    std::optional<double> max;
    double independent_value;
    double value_expected;
    double dValue_expected;
};

TEST(MaterialPropertyLib, Linear)
{
    double const y_ref = 1.0;
    double const m = 1.0;
    double const x_ref = 293.15;
    double const T = 303.15;

    std::vector<MaterialPropertyLibLinearTestParams> const params{
        {std::nullopt, std::nullopt, T, y_ref * (1 + m * (T - x_ref)),
         y_ref * m},
        // cutoffs but not effective
        {303., 304., T, y_ref * (1 + m * (T - x_ref)), y_ref * m},
        // min cutoff
        {304.15, std::nullopt, T, y_ref * (1 + m * (304.15 - x_ref)), 0.},
        // max cutoff
        {std::nullopt, 302.15, T, y_ref * (1 + m * (302.15 - x_ref)), 0.},
        // both cutoffs
        {304.15, 305.15, T, y_ref * (1 + m * (304.15 - x_ref)), 0.},
        {302.15, 302.15, T, y_ref * (1 + m * (302.15 - x_ref)), 0.}};

    auto value_or_none = [](std::optional<double> v) -> std::string
    { return v ? std::to_string(*v) : "(nullopt)"; };

    for (auto const& [min, max, T, value_expected, dValue_expected] : params)
    {
        INFO("T = {}, min = {}, max = {}", T, value_or_none(min),
             value_or_none(max));

        MaterialPropertyLib::IndependentVariable const iv{
            MaterialPropertyLib::Variable::temperature, x_ref, m, min, max};

        std::vector<MaterialPropertyLib::IndependentVariable> ivs{iv};
        MaterialPropertyLib::Linear linear_property{"linear", y_ref, ivs};

        MaterialPropertyLib::VariableArray variable_array;
        variable_array.temperature = T;
        ParameterLib::SpatialPosition const pos;
        double const time = std::numeric_limits<double>::quiet_NaN();
        double const dt = std::numeric_limits<double>::quiet_NaN();
        ASSERT_NEAR(std::get<double>(
                        linear_property.value(variable_array, pos, time, dt)),
                    value_expected,
                    1.e-10);
        ASSERT_EQ(std::get<double>(linear_property.dValue(
                      variable_array,
                      MaterialPropertyLib::Variable::liquid_phase_pressure, pos,
                      time, dt)),
                  0.0);
        ASSERT_NEAR(
            std::get<double>(linear_property.dValue(
                variable_array, MaterialPropertyLib::Variable::temperature, pos,
                time, dt)),
            dValue_expected, 1.e-16);
        ASSERT_EQ(
            std::get<double>(linear_property.d2Value(
                variable_array, MaterialPropertyLib::Variable::temperature,
                MaterialPropertyLib::Variable::temperature, pos, time, dt)),
            0.0);
    }
}

TEST(MaterialPropertyLib, LinearWithTime)
{
    double const y_ref = 1.0;
    double const m = 1.0;
    double const x_ref = 0.0;
    double const t = 10;

    std::vector<MaterialPropertyLibLinearTestParams> const params{
        {std::nullopt, std::nullopt, t, y_ref * (1 + m * (t - x_ref)),
         y_ref * m},
        // cutoffs but not effective
        {9, 11, t, y_ref * (1 + m * (t - x_ref)), y_ref * m},
        // min cutoff
        {11, std::nullopt, t, y_ref * (1 + m * (11 - x_ref)), 0.},
        // max cutoff
        {std::nullopt, 9, t, y_ref * (1 + m * (9 - x_ref)), 0.},
        // both cutoffs
        {11, 12, t, y_ref * (1 + m * (11 - x_ref)), 0.},
        {9, 9, t, y_ref * (1 + m * (9 - x_ref)), 0.}};

    auto value_or_none = [](std::optional<double> v) -> std::string
    { return v ? std::to_string(*v) : "(nullopt)"; };

    for (auto const& [min, max, t, value_expected, dValue_expected] : params)
    {
        INFO("t = {}, min = {}, max = {}", t, value_or_none(min),
             value_or_none(max));

        MaterialPropertyLib::IndependentVariable const iv{"t", x_ref, m, min,
                                                          max};

        std::vector<MaterialPropertyLib::IndependentVariable> ivs{iv};
        MaterialPropertyLib::Linear linear_property{"linear", y_ref, ivs};

        MaterialPropertyLib::VariableArray variable_array;

        ParameterLib::SpatialPosition const pos;
        double const dt = std::numeric_limits<double>::quiet_NaN();
        ASSERT_NEAR(
            std::get<double>(linear_property.value(variable_array, pos, t, dt)),
            value_expected,
            1.e-10);
        ASSERT_EQ(std::get<double>(linear_property.dValue(
                      variable_array,
                      MaterialPropertyLib::Variable::liquid_phase_pressure, pos,
                      t, dt)),
                  0.0);
    }
}

TEST(MaterialPropertyLib, LinearWithPos)
{
    double const y_ref = 1.0;
    double const m = 1.0;
    Eigen::Vector3d x_refs{50, 20, -10};
    double const coord = 10;

    MaterialPropertyLib::VariableArray const variable_array;

    double const time = 10;
    double const dt = std::numeric_limits<double>::quiet_NaN();

    std::vector<std::string> xyz{"x", "y", "z"};

    for (auto const& [component, dir] : ranges::views::enumerate(xyz))
    {
        INFO("component = {}, dir = {}", component, dir);

        double const x_ref = x_refs[component];

        std::vector<MaterialPropertyLibLinearTestParams> const params{
            {std::nullopt, std::nullopt, coord,
             y_ref * (1 + m * (coord - x_ref)), y_ref * m},
            // cutoffs but not effective
            {9, 11, coord, y_ref * (1 + m * (coord - x_ref)), y_ref * m},
            // min cutoff
            {11, std::nullopt, coord, y_ref * (1 + m * (11 - x_ref)), 0.},
            // max cutoff
            {std::nullopt, 9, coord, y_ref * (1 + m * (9 - x_ref)), 0.},
            // both cutoffs
            {11, 12, coord, y_ref * (1 + m * (11 - x_ref)), 0.},
            {9, 9, coord, y_ref * (1 + m * (9 - x_ref)), 0.}};

        auto value_or_none = [](std::optional<double> v) -> std::string
        { return v ? std::to_string(*v) : "(nullopt)"; };

        for (auto const& [min, max, x, value_expected, dValue_expected] :
             params)
        {
            INFO("coord = {}, min = {}, max = {}", x, value_or_none(min),
                 value_or_none(max));

            MaterialPropertyLib::IndependentVariable const iv{dir, x_ref, m,
                                                              min, max};

            std::vector<MaterialPropertyLib::IndependentVariable> ivs{iv};
            MaterialPropertyLib::Linear linear_property{"linear", y_ref, ivs};

            std::array<double, 3> coords = {-1, -2, -3};
            coords[component] = coord;
            ParameterLib::SpatialPosition const pos =
                ParameterLib::SpatialPosition{
                    {0}, {0}, MathLib::Point3d{coords}};

            EXPECT_NEAR(std::get<double>(linear_property.value(variable_array,
                                                               pos, time, dt)),
                        value_expected,
                        1.e-10);
            EXPECT_EQ(std::get<double>(linear_property.dValue(
                          variable_array,
                          MaterialPropertyLib::Variable::liquid_phase_pressure,
                          pos, time, dt)),
                      0.0);
        }
    }
}

TEST(MaterialPropertyLib, BiLinearTemperaturePosition)
{
    double const y_ref = 1.0;
    double const m = 1.0;
    double const x_ref = 0.0;
    double const T_ref = 273.15;
    double const T = 293.15;

    MaterialPropertyLib::IndependentVariable const iv_x{
        "x", x_ref, m, std::nullopt, std::nullopt};
    MaterialPropertyLib::IndependentVariable const iv_T{
        MaterialPropertyLib::Variable::temperature, T_ref, m, std::nullopt,
        std::nullopt};

    MaterialPropertyLib::Linear p{"bi-linear x T", y_ref, {iv_x, iv_T}};

    MaterialPropertyLib::VariableArray variable_array;
    variable_array.temperature = T;

    std::array<double, 3> const coords = {50, 20, -10};
    ParameterLib::SpatialPosition const pos =
        ParameterLib::SpatialPosition{{0}, {0}, MathLib::Point3d{coords}};
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    ASSERT_NEAR(std::get<double>(p.value(variable_array, pos, time, dt)),
                y_ref * (1. + m * (coords[0] - x_ref) + m * (T - T_ref)),
                1.e-10);
    ASSERT_EQ(std::get<double>(p.dValue(
                  variable_array, MaterialPropertyLib::Variable::temperature,
                  pos, time, dt)),
              y_ref * m);
    ASSERT_EQ(std::get<double>(
                  p.dValue(variable_array,
                           MaterialPropertyLib::Variable::liquid_phase_pressure,
                           pos, time, dt)),
              0.0);
}
