/*!
   \file
   \brief Test classes for water saturated enthalpy models.

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

#include "MaterialLib/MPL/Properties/Enthalpy/CreateWaterLiquidEnthalpyIAPWSIF97Region4.h"
#include "MaterialLib/MPL/Properties/Enthalpy/CreateWaterVapourEnthalpyIAPWSIF97Region4.h"
#include "TestMPL.h"

TEST(Material, checkWaterLiquidEnthalpyIAPWSIF97Region4)
{
    const char xml[] =
        "<property>"
        "   <name>saturation_enthalpy</name>"
        "   <type>WaterLiquidEnthalpyIAPWSIF97Region4</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml,
            MaterialPropertyLib::createWaterLiquidEnthalpyIAPWSIF97Region4);
    MaterialPropertyLib::Property const& property = *property_ptr;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    double const p[] = {611.213, 1.e6, 10.e6, 22.064e6};

    double const expected_h[] = {-41.555230442141529, 762647.04245918128,
                                 1407801.4124185159, 1923760.6606565041};

    for (int i = 0; i < 4; i++)
    {
        variable_array.liquid_phase_pressure = p[i];
        ASSERT_NEAR(expected_h[i],
                    property.template value<double>(variable_array, pos, t, dt),
                    1.e-6);
    }
}

TEST(Material, checkWaterVapourEnthalpyIAPWSIF97Region4)
{
    const char xml[] =
        "<property>"
        "   <name>saturation_enthalpy</name>"
        "   <type>WaterVapourEnthalpyIAPWSIF97Region4</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml,
            MaterialPropertyLib::createWaterVapourEnthalpyIAPWSIF97Region4);
    MaterialPropertyLib::Property const& property = *property_ptr;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();
    double const p[] = {611.213, 1.e6, 10.e6, 22.064e6};

    double const expected_h[] = {2500.775234200720, 2776.989174073231,
                                 2725.344627240651, 2303.146272347803};

    for (int i = 0; i < 4; i++)
    {
        variable_array.liquid_phase_pressure = p[i];
        ASSERT_NEAR(
            expected_h[i],
            property.template value<double>(variable_array, pos, t, dt) / 1e3,
            1.e-8);
    }
}
