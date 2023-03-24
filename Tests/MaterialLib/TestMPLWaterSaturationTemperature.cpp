/*!
   \file
   \brief Test classes for fluid density models.

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

#include "MaterialLib/MPL/Properties/CreateWaterSaturationTemperatureIAPWSIF97Region4.h"
#include "TestMPL.h"

TEST(Material, checkSaturationTemperatureIAPWSIF97Region4)
{
    const char xml[] =
        "<property>"
        "   <name>Saturation_Temperature</name>"
        "   <type>WaterSaturationTemperatureIAPWSIF97Region4</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::
                     createWaterSaturationTemperatureIAPWSIF97Region4);
    MaterialPropertyLib::Property const& property = *property_ptr;

    const double p = 5.e+6;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    variable_array.phase_pressure = p;

    const double T_s = 537.09287118633063;

    ASSERT_NEAR(T_s,
                property.template value<double>(variable_array, pos, t, dt),
                1.e-10);
}
