// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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

    double const expected_h[] = {-41.55718122033141, 762682.8443354114,
                                 1407867.500568234, 1923850.9701145296};

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

    double const expected_h[] = {2500.8926311621317, 2777.1195376846626,
                                 2725.4725664387292, 2303.2543917702774};

    for (int i = 0; i < 4; i++)
    {
        variable_array.liquid_phase_pressure = p[i];
        ASSERT_NEAR(
            expected_h[i],
            property.template value<double>(variable_array, pos, t, dt) / 1e3,
            1.e-8);
    }
}
