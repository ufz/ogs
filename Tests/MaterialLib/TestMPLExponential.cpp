// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <autocheck/autocheck.hpp>
#include <cmath>

#include "MaterialLib/MPL/Properties/Exponential.h"
#include "Tests/AutoCheckTools.h"

namespace ac = autocheck;
namespace MPL = MaterialPropertyLib;

struct MaterialPropertyLibExponentialProperty : public ::testing::Test
{
    void SetUp() override
    {
        double const y_offset = 1.;
        double const y_ref = 1.;
        double const T_ref = 1.;
        double const factor = 2.;
        p = std::make_unique<MPL::Exponential>(
            "exponential", y_offset, y_ref,
            MPL::ExponentData{MPL::Variable::temperature, T_ref, factor});
    }
    std::unique_ptr<MPL::Property> p;
    ac::gtest_reporter gtest_reporter;
};

// First order derivative approximated with second order central differences.
template <typename F>
double dydx_C2(double const x, F&& f, double const h)
{
    return (f(x + h) - f(x - h)) / (2 * h);
}
// Second order derivative approximated with second order central differences.
template <typename F>
double d2ydx2_C2(double const x, F&& f, double const h)
{
    return (f(x + h) - 2 * f(x) + f(x - h)) / (h * h);
}

TEST_F(MaterialPropertyLibExponentialProperty, TestNumericalDerivatives)
{
    double const eps = 1e-5;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    auto y = [&](double const T)
    {
        MPL::VariableArray variable_array;
        variable_array.temperature = T;
        return p->template value<double>(variable_array, pos, time, dt);
    };

    auto f = [&](double const T)
    {
        MPL::VariableArray variable_array;
        variable_array.temperature = T;
        double const v =
            p->template value<double>(variable_array, pos, time, dt);
        double const dv = p->template dValue<double>(
            variable_array, MPL::Variable::temperature, pos, time, dt);
        double const dv2 = p->template d2Value<double>(
            variable_array, MPL::Variable::temperature,
            MPL::Variable::temperature, pos, time, dt);

        double const Dv = dydx_C2(T, y, eps);
        double const Dv2 = d2ydx2_C2(T, y, eps);

        if ((std::abs(dv - Dv) > 1e-9 * v) ||
            (std::abs(dv2 - Dv2) > 1.5e-4 * v))
            INFO("{} {} {}", T, std::abs(dv - Dv) / v, std::abs(dv2 - Dv2) / v);

        return (std::abs(dv - Dv) <= 1e-9 * v) &&
               (std::abs(dv2 - Dv2) <= 1.5e-4 * v);
    };

    // Limit values to avoid +-inf.
    auto gen = ac::IntervalGenerator(-100., 100.);
    ac::check<double>(f, 10000, ac::make_arbitrary(gen), gtest_reporter);
}

TEST(MaterialPropertyLib, Exponential)
{
    double const y_offset = -1e-3;
    double const y_ref = 1e-3;
    double const reference_condition = 20.0;
    double const factor = 1 / 75.0;
    MPL::ExponentData const exp_data{MPL::Variable::temperature,
                                     reference_condition, factor};
    MPL::Property const& p =
        MPL::Exponential{"exponential", y_offset, y_ref, exp_data};

    double const T = 20.;
    MPL::VariableArray variable_array;
    variable_array.temperature = T;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    ASSERT_NEAR(
        p.template value<double>(variable_array, pos, time, dt),
        y_offset + y_ref * (std::exp(factor * (T - reference_condition))),
        1.e-10);
    ASSERT_EQ(p.template dValue<double>(variable_array,
                                        MPL::Variable::liquid_phase_pressure,
                                        pos, time, dt),
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

TEST(MaterialPropertyLib, ExponentialTime)
{
    double const y_offset = -1e-3;
    double const y_ref = 1e-3;
    double const reference_condition = 0.0;
    double const factor = 1 / 75.0;
    MPL::ExponentData const exp_data{"t", reference_condition, factor};
    MPL::Property const& p =
        MPL::Exponential{"exponential", y_offset, y_ref, exp_data};

    MPL::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const time = 20;
    double const dt = std::numeric_limits<double>::quiet_NaN();
    ASSERT_NEAR(
        p.template value<double>(variable_array, pos, time, dt),
        y_offset + y_ref * (std::exp(factor * (time - reference_condition))),
        1.e-10);
    ASSERT_EQ(p.template dValue<double>(variable_array,
                                        MPL::Variable::liquid_phase_pressure,
                                        pos, time, dt),
              0.0);
    ASSERT_EQ(p.template d2Value<double>(
                  variable_array, MPL::Variable::liquid_phase_pressure,
                  MPL::Variable::liquid_phase_pressure, pos, time, dt),
              0.0);
}

TEST(MaterialPropertyLib, ExponentialPos)
{
    double const y_offset = -1e-3;
    double const y_ref = 1e-3;
    double const reference_condition = 0.0;
    double const factor = 1 / 75.0;
    MPL::ExponentData const exp_data_x{"x", reference_condition, factor};
    MPL::ExponentData const exp_data_y{"y", reference_condition, factor};
    MPL::ExponentData const exp_data_z{"z", reference_condition, factor};

    MPL::Property const& p_x =
        MPL::Exponential{"exponential", y_offset, y_ref, exp_data_x};
    MPL::Property const& p_y =
        MPL::Exponential{"exponential", y_offset, y_ref, exp_data_y};
    MPL::Property const& p_z =
        MPL::Exponential{"exponential", y_offset, y_ref, exp_data_z};

    MPL::VariableArray variable_array;
    std::array<double, 3> coords = {0, 20, 0};
    ParameterLib::SpatialPosition const pos =
        ParameterLib::SpatialPosition{{0}, {0}, MathLib::Point3d{coords}};
    double const time = 20;
    double const dt = std::numeric_limits<double>::quiet_NaN();

    ASSERT_NEAR(
        p_x.template value<double>(variable_array, pos, time, dt),
        y_offset +
            y_ref * (std::exp(factor * (coords[0] - reference_condition))),
        1.e-10);
    ASSERT_EQ(p_x.template dValue<double>(variable_array,
                                          MPL::Variable::liquid_phase_pressure,
                                          pos, time, dt),
              0.0);
    ASSERT_EQ(p_x.template d2Value<double>(
                  variable_array, MPL::Variable::liquid_phase_pressure,
                  MPL::Variable::liquid_phase_pressure, pos, time, dt),
              0.0);

    ASSERT_NEAR(
        p_y.template value<double>(variable_array, pos, time, dt),
        y_offset +
            y_ref * (std::exp(factor * (coords[1] - reference_condition))),
        1.e-10);
    ASSERT_EQ(p_y.template dValue<double>(variable_array,
                                          MPL::Variable::liquid_phase_pressure,
                                          pos, time, dt),
              0.0);
    ASSERT_EQ(p_y.template d2Value<double>(
                  variable_array, MPL::Variable::liquid_phase_pressure,
                  MPL::Variable::liquid_phase_pressure, pos, time, dt),
              0.0);

    ASSERT_NEAR(
        p_z.template value<double>(variable_array, pos, time, dt),
        y_offset +
            y_ref * (std::exp(factor * (coords[2] - reference_condition))),
        1.e-10);
    ASSERT_EQ(p_z.template dValue<double>(variable_array,
                                          MPL::Variable::liquid_phase_pressure,
                                          pos, time, dt),
              0.0);
    ASSERT_EQ(p_z.template d2Value<double>(
                  variable_array, MPL::Variable::liquid_phase_pressure,
                  MPL::Variable::liquid_phase_pressure, pos, time, dt),
              0.0);
}
