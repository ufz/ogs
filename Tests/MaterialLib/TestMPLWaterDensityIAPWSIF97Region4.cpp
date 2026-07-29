// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <memory>

#include "MaterialLib/MPL/Properties/Density/CreateWaterLiquidDensityIAPWSIF97Region4.h"
#include "MaterialLib/MPL/Properties/Density/CreateWaterVapourDensityIAPWSIF97Region4.h"
#include "TestMPL.h"

TEST(Material, checkWaterLiquidDensityIAPWSIF97Region4)
{
    const char xml[] =
        "<property>"
        "   <name>saturation_density</name>"
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

    double const p[] = {611.213, 1.e6, 10.e6, 22.064e6};

    double const expected_rho[] = {999.7930660007243, 887.1274516747791,
                                   688.4113330921692, 445.206163508301};

    for (int i = 0; i < 4; i++)
    {
        variable_array.liquid_phase_pressure = p[i];
        ASSERT_NEAR(expected_rho[i],
                    property.template value<double>(variable_array, pos, t, dt),
                    1.e-9);
    }
}

TEST(Material, checkWaterVapourDensityIAPWSIF97Region4)
{
    const char xml[] =
        "<property>"
        "   <name>saturation_density</name>"
        "   <type>WaterVapourDensityIAPWSIF97Region4</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createWaterVapourDensityIAPWSIF97Region4);
    MaterialPropertyLib::Property const& property = *property_ptr;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    double const p[] = {611.213, 1.e6, 10.e6, 22.064e6};

    double const expected_rho[] = {0.004851081195212619, 5.14538585318268,
                                   55.45212134316541, 224.57431856923355};

    for (int i = 0; i < 4; i++)
    {
        variable_array.liquid_phase_pressure = p[i];
        ASSERT_NEAR(expected_rho[i],
                    property.template value<double>(variable_array, pos, t, dt),
                    1.e-9);
    }
}
