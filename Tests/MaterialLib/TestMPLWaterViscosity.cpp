/**
 *  \file
 *  \brief Test viscosity models
 *
 *  \copyright
 *   Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */

#include <gtest/gtest.h>

#include <memory>

#include "MaterialLib/MPL/Properties/Viscosity/CreateLiquidViscosityVogels.h"
#include "MaterialLib/MPL/Properties/Viscosity/CreateWaterViscosityIAPWS.h"
#include "TestMPL.h"

TEST(Material, checkWaterViscosityIAPWS_)
{
    const char xml[] =
        "<property>"
        "  <name>viscosity</name>"
        "  <type>WaterViscosityIAPWS</type>"
        "</property>";

    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr =
        Tests::createTestProperty(
            xml, MaterialPropertyLib::createWaterViscosityIAPWS);
    MaterialPropertyLib::Property const& property = *property_ptr;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    // Test data provided on http://www.iapws.org/relguide/visc.pdf
    const double T[] = {298.15, 298.15, 373.15,  433.15,  433.15, 873.15,
                        873.15, 873.15, 1173.15, 1173.15, 1173.15};
    const double rho[] = {998, 1200, 1000, 1, 1000, 1, 100, 600, 1, 100, 400};

    const double expected_mumu[] = {
        889.735100, 1437.649467, 307.883622, 14.538324, 217.685358, 32.619287,
        35.802262,  77.430195,   44.217245,  47.640433, 64.154608};

    const double perturbation = 1.e-9;
    for (int i = 0; i < 11; i++)
    {
        // Test mu
        variable_array.temperature = T[i];
        variable_array.density = rho[i];
        const double mu =
            property.template value<double>(variable_array, pos, t, dt);
        ASSERT_NEAR(expected_mumu[i] * 1.e-6, mu, 1.e-9);

        const double dmu_dT = property.template dValue<double>(
            variable_array, MaterialPropertyLib::Variable::temperature, pos, t,
            dt);
        const double dmu_drho = property.template dValue<double>(
            variable_array, MaterialPropertyLib::Variable::density, pos, t, dt);

        // Test dmu/dT
        variable_array.temperature = T[i] + perturbation;
        double mu1 =
            property.template value<double>(variable_array, pos, t, dt);
        ASSERT_NEAR((mu1 - mu) / perturbation, dmu_dT, 8.e-6);

        // Test dmu/drho
        variable_array.temperature = T[i];
        variable_array.density = rho[i] + perturbation;
        mu1 = property.template value<double>(variable_array, pos, t, dt);

        ASSERT_NEAR((mu1 - mu) / perturbation, dmu_drho, 8.e-6);
    }
}
TEST(Material, checkVogelViscosity_)
{
    const char xml_w[] =
        "<property>"
        "  <name>viscosity</name>"
        "  <type>LiquidViscosityVogels</type>"
        "  <liquid_type>Water</liquid_type>"
        "</property>";
    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr_w =
        Tests::createTestProperty(
            xml_w, MaterialPropertyLib::createLiquidViscosityVogels);
    MaterialPropertyLib::Property const& property_w = *property_ptr_w;

    MaterialPropertyLib::VariableArray variable_array;
    ParameterLib::SpatialPosition const pos;
    double const t = std::numeric_limits<double>::quiet_NaN();
    double const dt = std::numeric_limits<double>::quiet_NaN();

    variable_array.temperature = 303.;
    auto const mu_w =
        property_w.template value<double>(variable_array, pos, t, dt);
    auto const dmu_w_dT = property_w.template dValue<double>(
        variable_array, MaterialPropertyLib::Variable::temperature, pos, t, dt);
    ASSERT_NEAR(0.802657e-3, mu_w, 1.e-5);
    ASSERT_NEAR(-1.87823e-5, dmu_w_dT, 1.e-5);

    const char xml_co2[] =
        "<property>"
        "  <name>viscosity</name>"
        "  <type>LiquidViscosityVogels</type>"
        "  <liquid_type> CO2 </liquid_type>"
        "</property>";
    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr_co2 =
        Tests::createTestProperty(
            xml_co2, MaterialPropertyLib::createLiquidViscosityVogels);
    MaterialPropertyLib::Property const& property_co2 = *property_ptr_co2;

    variable_array.temperature = 255.04;
    auto const mu_co2 =
        property_co2.template value<double>(variable_array, pos, t, dt);
    auto const dmu_co2_dT = property_co2.template dValue<double>(
        variable_array, MaterialPropertyLib::Variable::temperature, pos, t, dt);
    ASSERT_NEAR(0.137956e-3, mu_co2, 1.e-5);
    ASSERT_NEAR(-2.35664e-6, dmu_co2_dT, 1.e-5);

    const char xml_ch4[] =
        "<property>"
        "  <name>viscosity</name>"
        "  <type>LiquidViscosityVogels</type>"
        "  <liquid_type> CH4 </liquid_type>"
        "</property>";
    std::unique_ptr<MaterialPropertyLib::Property> const property_ptr_ch4 =
        Tests::createTestProperty(
            xml_ch4, MaterialPropertyLib::createLiquidViscosityVogels);
    MaterialPropertyLib::Property const& property_ch4 = *property_ptr_ch4;

    variable_array.temperature = 172.0;
    auto const mu_ch4 =
        property_ch4.template value<double>(variable_array, pos, t, dt);
    auto const dmu_ch4_dT = property_ch4.template dValue<double>(
        variable_array, MaterialPropertyLib::Variable::temperature, pos, t, dt);
    ASSERT_NEAR(0.352072e-4, mu_ch4, 1.e-5);
    ASSERT_NEAR(-2.35664e-6, dmu_ch4_dT, 1.e-5);
}
