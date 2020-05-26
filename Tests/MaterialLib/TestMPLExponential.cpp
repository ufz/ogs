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
#include <cmath>

#include <gtest/gtest.h>

#include "MaterialLib/MPL/Properties/Exponential.h"

TEST(MaterialPropertyLib, Exponential)
{
    double const y_ref = 1e-3;
    double const reference_condition = 20.0;
    double const factor = 1 / 75.0;
    MaterialPropertyLib::ExponentData const exp_data{
        MaterialPropertyLib::Variable::temperature, reference_condition, factor};
    MaterialPropertyLib::Exponential exp_property{"exponential", y_ref,
                                                  exp_data};

    MaterialPropertyLib::VariableArray variable_array;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::temperature)] = 20.;
    ParameterLib::SpatialPosition const pos;
    double const time = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    ASSERT_NEAR(
        std::get<double>(exp_property.value(variable_array, pos, time, dt)),
        y_ref * (std::exp(-factor *
                          (std::get<double>(variable_array[static_cast<int>(
                               MaterialPropertyLib::Variable::temperature)]) -
                           reference_condition))),
        1.e-10);
    ASSERT_EQ(std::get<double>(exp_property.dValue(
                  variable_array, MaterialPropertyLib::Variable::phase_pressure,
                  pos, time, dt)),
              0.0);
    ASSERT_NEAR(
        std::get<double>(exp_property.dValue(
            variable_array, MaterialPropertyLib::Variable::temperature, pos,
            time, dt)),
        -y_ref * factor *
            std::exp(-factor *
                     (std::get<double>(variable_array[static_cast<int>(
                          MaterialPropertyLib::Variable::temperature)]) -
                      reference_condition)),
        1.e-16);
    ASSERT_NEAR(
        std::get<double>(exp_property.d2Value(
            variable_array, MaterialPropertyLib::Variable::temperature,
            MaterialPropertyLib::Variable::temperature, pos, time, dt)),
        y_ref * std::pow(factor, 2) *
            std::exp(-factor *
                     (std::get<double>(variable_array[static_cast<int>(
                          MaterialPropertyLib::Variable::temperature)]) -
                      reference_condition)),
        1.e-16);
}

