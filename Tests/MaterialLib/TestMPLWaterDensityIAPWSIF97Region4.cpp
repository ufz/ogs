/*!
   \file
   \brief Test classes for water saturated density models.

   \author Chaofan Chen
   \date March 2023

   \copyright
    Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#include <gtest/gtest.h>

#include <memory>

#include "MaterialLib/MPL/Properties/Density/CreateWaterLiquidDensityIAPWSIF97Region4.h"
#include "TestMPL.h"

TEST(Material, checkWaterLiquidDensityIAPWSIF97Region4)
{
    const char xml[] =
        "<property>"
        "   <name>density</name>"
        "   <type>WaterLiquidDensityIAPWSIF97Region4</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createWaterLiquidDensityIAPWSIF97Region4);
    MaterialPropertyLib::Property const& property = *property_ptr;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    double const p[] = {611.213, 1.e+6, 10.e6};

    double const expected_rho[] = {999.84000051382395, 887.16909718765328,
                                   688.44365003037342};

    for (int i = 0; i < 3; i++)
    {
        variable_array.phase_pressure = p[i];
        ASSERT_NEAR(expected_rho[i],
                    property.template value<double>(variable_array, pos, t, dt),
                    1.e-12);
    }
}
