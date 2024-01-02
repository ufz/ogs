/*!
   \file
   \brief Test class for water temperature in IAPWSIF97 region1 model.

   \author Chaofan Chen
   \date March 2023

   \copyright
    Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#include <gtest/gtest.h>

#include <memory>

#include "MaterialLib/MPL/Properties/CreateWaterTemperatureIAPWSIF97Region1.h"
#include "TestMPL.h"

TEST(Material, checkWaterTemperatureIAPWSIF97Region1)
{
    const char xml[] =
        "<property>"
        "   <name>temperature</name>"
        "   <type>WaterTemperatureIAPWSIF97Region1</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createWaterTemperatureIAPWSIF97Region1);
    MaterialPropertyLib::Property const& property = *property_ptr;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    double const p[] = {0, 3.e6, 80.e6, 80.e6};
    double const h[] = {0, 5.e5, 5.e5, 1.5e6};

    double const expected_T[] = {0.273138540e3, 0.391798509e3, 0.378108626e3,
                                 0.611041229e3};

    for (int i = 0; i < 4; i++)
    {
        variable_array.liquid_phase_pressure = p[i];
        variable_array.enthalpy = h[i];

        ASSERT_NEAR(expected_T[i],
                    property.template value<double>(variable_array, pos, t, dt),
                    1.e-6);
    }
}
