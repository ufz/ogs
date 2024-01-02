/*!
   \file
   \brief Test classes for fluid density models.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#include <gtest/gtest.h>

#include <memory>

#include "MaterialLib/MPL/Properties/Density/CreateWaterDensityIAPWSIF97Region1.h"
#include "TestMPL.h"

TEST(Material, checkWaterDensityIAPWSIF97Region1_)
{
    const char xml[] =
        "<property>"
        "   <name>density</name>"
        "   <type>WaterDensityIAPWSIF97Region1</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createWaterDensityIAPWSIF97Region1);
    MaterialPropertyLib::Property const& property = *property_ptr;

    const double T = 473.15;
    const double p = 4.e+7;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    variable_array.temperature = T;
    variable_array.liquid_phase_pressure = p;

    const double rho_expected = 890.98496087498256;

    ASSERT_NEAR(rho_expected,
                property.template value<double>(variable_array, pos, t, dt),
                1.e-10);

    const double drho_dT = property.template dValue<double>(
        variable_array, MaterialPropertyLib::Variable::temperature, pos, t, dt);
    const double drho_dp = property.template dValue<double>(
        variable_array, MaterialPropertyLib::Variable::liquid_phase_pressure,
        pos, t, dt);

    const double perturbation = 1.e-4;

    // Test the differentiation: with respect to temperature:
    variable_array.temperature += perturbation;
    const double rho_T1 =
        property.template value<double>(variable_array, pos, t, dt);
    ASSERT_NEAR((rho_T1 - rho_expected) / perturbation, drho_dT, 1.e-6);

    // Test the differentiation: with respect to pressure:
    variable_array.liquid_phase_pressure += perturbation;
    const double rho_p1 =
        property.template value<double>(variable_array, pos, t, dt);
    ASSERT_NEAR((rho_p1 - rho_T1) / perturbation, drho_dp, 1.e-6);
}
