/**
 *  \brief Test viscosity models
 *
 *  \copyright
 *   Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file
 *
 */

#include <gtest/gtest.h>

#include <memory>
#include <cmath>

#include "Tests/TestTools.h"

#include "MaterialLib/Fluid/Viscosity/CreateViscosityModel.h"
#include "MaterialLib/PhysicalConstant.h"

using namespace MaterialLib;
using namespace MaterialLib::Fluid;
using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

std::unique_ptr<FluidProperty> createTestViscosityModel(const char xml[])
{
    auto const ptree = Tests::readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("viscosity");
    return MaterialLib::Fluid::createViscosityModel(sub_config);
}

TEST(Material, checkConstantViscosity)
{
    const char xml[] =
        "<viscosity>"
        "   <type>Constant</type>"
        "   <value> 1.e-4 </value> "
        "</viscosity>";
    auto const mu = createTestViscosityModel(xml);

    ArrayType dummy;
    ASSERT_EQ(1.e-4, mu->getValue(dummy));
    ASSERT_EQ(
        0.0, mu->getdValue(dummy, MaterialLib::Fluid::PropertyVariableType::T));
}

TEST(Material, checkTemperatureDependentViscosity)
{
    const char xml[] =
        "<viscosity>"
        "  <type>TemperatureDependent</type>"
        "  <mu0>1.e-3 </mu0>"
        "   <tc>293.</tc>"
        "   <tv>368.</tv>"
        "</viscosity>";
    auto const mu = createTestViscosityModel(xml);

    ArrayType vars;
    vars[0] = 350.0;
    const double mu_expected = 1.e-3 * std::exp(-(vars[0] - 293) / 368);
    ASSERT_NEAR(mu_expected, mu->getValue(vars), 1.e-10);
    const double dmu_dT_expected =
        -1.e-3/368 * std::exp(-(vars[0] - 293) / 368);
    ASSERT_NEAR(
        dmu_dT_expected,
        mu->getdValue(vars, MaterialLib::Fluid::PropertyVariableType::T),
        1.e-10);
}

TEST(Material, checkLinearPressureDependentViscosity)
{
    const char xml[] =
        "<viscosity>"
        "  <type>LinearPressure</type>"
        "  <mu0>1.e-3 </mu0>"
        "   <p0>1.e+5</p0>"
        "   <gamma>1.e-6</gamma>"
        "</viscosity>";
    const auto mu = createTestViscosityModel(xml);

    ArrayType vars;
    vars[0] = 293.;
    vars[1] = 2.e+6;
    ASSERT_NEAR(1.e-3 * (1. + 1.e-6 * (vars[1] - 1.e+5)), mu->getValue(vars),
                1.e-10);
    ASSERT_NEAR(
        1.e-9,
        mu->getdValue(vars, MaterialLib::Fluid::PropertyVariableType::p),
        1.e-10);
}

TEST(Material, checkVogelViscosity)
{
    const char xml_w[] =
        "<viscosity>"
        "  <type>Vogels</type>"
        "  <liquid_type>Water </liquid_type>"
        "</viscosity>";
    const auto mu_w = createTestViscosityModel(xml_w);
    ArrayType vars;
    vars[0] = 303.;
    const auto var_type = MaterialLib::Fluid::PropertyVariableType::T;
    ASSERT_NEAR(0.802657e-3, mu_w->getValue(vars), 1.e-5);
    ASSERT_NEAR(-1.87823e-5, mu_w->getdValue(vars, var_type), 1.e-5);

    const char xml_co2[] =
        "<viscosity>"
        "  <type>Vogels</type>"
        "  <liquid_type> CO2 </liquid_type>"
        "</viscosity>";
    const auto mu_co2 = createTestViscosityModel(xml_co2);
    vars[0] = 255.04;
    ASSERT_NEAR(0.137956e-3, mu_co2->getValue(vars), 1.e-5);
    ASSERT_NEAR(-2.35664e-6, mu_co2->getdValue(vars, var_type), 1.e-5);

    const char xml_ch4[] =
        "<viscosity>"
        "  <type>Vogels</type>"
        "  <liquid_type> CH4 </liquid_type>"
        "</viscosity>";
    const auto mu_ch4 = createTestViscosityModel(xml_ch4);
    vars[0] = 172.0;
    ASSERT_NEAR(0.352072e-4, mu_ch4->getValue(vars), 1.e-5);
    ASSERT_NEAR(-2.35664e-6, mu_ch4->getdValue(vars, var_type), 1.e-5);
}

TEST(Material, checkWaterViscosityIAPWS)
{
    const char xml_w[] =
        "<viscosity>"
        "  <type>WaterViscosityIAPWS</type>"
        "</viscosity>";

    // Test data provided on http://www.iapws.org/relguide/visc.pdf
    const auto mu_w = createTestViscosityModel(xml_w);
    const double T[] = {298.15, 298.15, 373.15,  433.15,  433.15, 873.15,
                        873.15, 873.15, 1173.15, 1173.15, 1173.15};
    const double rho[] = {998, 1200, 1000, 1, 1000, 1, 100, 600, 1, 100, 400};

    const double expected_mumu[] = {
        889.735100, 1437.649467, 307.883622, 14.538324, 217.685358, 32.619287,
        35.802262,  77.430195,   44.217245,  47.640433, 64.154608};

    const double perturbation = 1.e-9;
    ArrayType vars;
    for (int i = 0; i < 11; i++)
    {
        // Test mu
        vars[static_cast<unsigned>(PropertyVariableType::T)] = T[i];
        vars[static_cast<unsigned>(PropertyVariableType::rho)] = rho[i];
        const double mu = mu_w->getValue(vars);
        ASSERT_NEAR(expected_mumu[i] * 1.e-6, mu, 1.e-9);

        const double dmu_dT = mu_w->getdValue(vars, PropertyVariableType::T);
        const double dmu_drho =
            mu_w->getdValue(vars, PropertyVariableType::rho);

        // Test dmu/dT
        vars[static_cast<unsigned>(PropertyVariableType::T)] =
            T[i] + perturbation;
        double mu1 = mu_w->getValue(vars);
        ASSERT_NEAR((mu1 - mu) / perturbation, dmu_dT, 1.e-7);

        // Test dmu/drho
        vars[static_cast<unsigned>(PropertyVariableType::T)] = T[i];
        vars[static_cast<unsigned>(PropertyVariableType::rho)] =
            rho[i] + perturbation;
        mu1 = mu_w->getValue(vars);

        ASSERT_NEAR((mu1 - mu) / perturbation, dmu_drho, 1.e-7);
    }
}
