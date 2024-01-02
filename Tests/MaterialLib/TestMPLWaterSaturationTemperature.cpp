/*!
   \file
   \brief Test classes for fluid density models.

   \author Chaofan Chen
   \date March 2023

   \copyright
    Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#include <gtest/gtest.h>

#include <array>
#include <memory>

#include "MaterialLib/MPL/Properties/CreateWaterSaturationTemperatureIAPWSIF97Region4.h"
#include "TestMPL.h"

TEST(Material, checkSaturationTemperatureIAPWSIF97Region4)
{
    const char xml[] =
        "<property>"
        "   <name>saturation_temperature</name>"
        "   <type>WaterSaturationTemperatureIAPWSIF97Region4</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::
                     createWaterSaturationTemperatureIAPWSIF97Region4);
    MaterialPropertyLib::Property const& property = *property_ptr;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    std::array const pressures = {611.213, 0.1e6, 1.e6, 10.e6, 22.064e6};

    std::array const expected_Ts = {273.150007261, 0.372755919e3, 0.453035632e3,
                                    0.584149488e3, 647.095999998};

    for (std::size_t i = 0; i < pressures.size(); i++)
    {
        variable_array.liquid_phase_pressure = pressures[i];
        ASSERT_NEAR(expected_Ts[i],
                    property.template value<double>(variable_array, pos, t, dt),
                    1.e-6);
    }
}
