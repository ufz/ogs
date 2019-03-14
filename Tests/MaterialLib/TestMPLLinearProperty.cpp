/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <gtest/gtest.h>

#include "MaterialLib/MPL/Properties/LinearProperty.h"

TEST(MaterialPropertyLib, LinearProperty)
{
    double const y_ref = 1.0;
    double const m = 1.0;
    double const x_ref = 293.15;
    MaterialPropertyLib::IndependentVariable const iv{
        MaterialPropertyLib::Variable::temperature, x_ref, m};
    MaterialPropertyLib::LinearProperty linear_property{y_ref, iv};

    MaterialPropertyLib::VariableArray variable_array;
    variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::temperature)] = 303.15;
    ASSERT_NEAR(
        boost::get<double>(linear_property.value(variable_array)),
        y_ref * (1 + m * (boost::get<double>(variable_array[static_cast<int>(
                              MaterialPropertyLib::Variable::temperature)]) -
                          x_ref)),
        1.e-10);
    ASSERT_EQ(
        boost::get<double>(linear_property.dValue(
            variable_array, MaterialPropertyLib::Variable::phase_pressure)),
        0.0);
    ASSERT_NEAR(
        boost::get<double>(linear_property.dValue(
            variable_array, MaterialPropertyLib::Variable::temperature)),
        y_ref * m, 1.e-16);
    ASSERT_EQ(boost::get<double>(linear_property.d2Value(
                  variable_array,
                  MaterialPropertyLib::Variable::temperature,
                  MaterialPropertyLib::Variable::temperature)),
              0.0);
}

