/*!
   \file  TestFluidDensity.cpp
   \brief Test classes for fluid density models.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#include <gtest/gtest.h>

#include <memory>

#include "Tests/TestTools.h"

#include "MaterialLib/Fluid/Density/CreateFluidDensityModel.h"
#include "MaterialLib/PhysicalConstant.h"

using namespace MaterialLib;
using namespace MaterialLib::Fluid;

using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;
//----------------------------------------------------------------------------
// Test density models.
std::unique_ptr<FluidProperty> createTestFluidDensityModel(const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("density");
    return MaterialLib::Fluid::createFluidDensityModel(sub_config);
}

TEST(Material, checkConstantDensity)
{
    const char xml[] =
        "<density>"
        "   <type>Constant</type>"
        "   <value> 998.0 </value> "
        "</density>";
    const auto rho = createTestFluidDensityModel(xml);

    ArrayType dummy;
    ASSERT_EQ(998.0, rho->getValue(dummy));
    ASSERT_EQ(
        0.0,
        rho->getdValue(dummy, MaterialLib::Fluid::PropertyVariableType::T));
}

TEST(Material, checkIdealGasLaw)
{
    const char xml[] =
        "<density>"
        "   <type>IdealGasLaw</type>"
        "   <molar_mass> 28.96 </molar_mass> "
        "</density>";
    const auto rho = createTestFluidDensityModel(xml);

    const double molar_air = 28.96;
    const double T = 290.;
    const double p = 1.e+5;
    const double R = PhysicalConstant::IdealGasConstant;
    const double expected_air_dens = molar_air * p / (R * T);
    ArrayType vars = {{290, 1.e+5}};
    ASSERT_NEAR(expected_air_dens, rho->getValue(vars), 1.e-10);

    const double expected_d_air_dens_dT = -molar_air * p / (R * T * T);
    ASSERT_NEAR(expected_d_air_dens_dT,
                rho->getdValue(vars, Fluid::PropertyVariableType::T), 1.e-10);

    const double expected_d_air_dens_dp = molar_air / (R * T);
    ASSERT_NEAR(expected_d_air_dens_dp,
                rho->getdValue(vars, Fluid::PropertyVariableType::p), 1.e-10);
}

TEST(Material, checkLinearTemperatureDependentDensity)
{
    const char xml[] =
        "<density>"
        "   <type>TemperatureDependent</type>"
        "   <temperature0> 293.0 </temperature0> "
        "   <beta> 4.3e-4 </beta> "
        "   <rho0>1000.</rho0>"
        "</density>";

    const auto rho = createTestFluidDensityModel(xml);

    ArrayType vars;
    vars[0] = 273.1;
    ASSERT_NEAR(1000.0 * (1 - 4.3e-4 * (vars[0] - 293.0)), rho->getValue(vars),
                1.e-10);
    ASSERT_NEAR(-1000.0 * 4.3e-4,
                rho->getdValue(vars, Fluid::PropertyVariableType::T), 1.e-10);
}

TEST(Material, checkLiquidDensity)
{
    const char xml[] =
        "<density>"
        "   <type>LiquidDensity</type>"
        "   <temperature0> 273.15 </temperature0> "
        "   <p0> 1.e+5 </p0> "
        "   <bulk_modulus> 2.15e+9 </bulk_modulus> "
        "   <beta> 2.0e-4 </beta> "
        "   <rho0>999.8</rho0>"
        "</density>";
    const auto rho = createTestFluidDensityModel(xml);

    const ArrayType vars = {{273.15 + 60.0, 1.e+6}};
    const double T0 = 273.15;
    const double p0 = 1.e+5;
    const double rho0 = 999.8;
    const double K = 2.15e+9;
    const double beta = 2.e-4;
    const double T = vars[0];
    const double p = vars[1];

    const double fac_T = 1. + beta * (T - T0);
    ASSERT_NEAR(rho0 / fac_T / (1. - (p - p0) / K), rho->getValue(vars),
                1.e-10);

    // Test the derivative with respect to temperature.
    ASSERT_NEAR(-beta * rho0 / (fac_T * fac_T) / (1. - (p - p0) / K),
                rho->getdValue(vars, Fluid::PropertyVariableType::T), 1.e-10);

    // Test the derivative with respect to pressure.
    const double fac_p = 1. - (p - p0) / K;
    ASSERT_NEAR(rho0 / (1. + beta * (T - T0)) / (fac_p * fac_p * K),
                rho->getdValue(vars, Fluid::PropertyVariableType::p), 1.e-10);
}

TEST(Material, checkWaterDensityIAPWSIF97Region1)
{
    const char xml[] =
        "<density>"
        "   <type>WaterDensityIAPWSIF97Region1</type>"
        "</density>";
    const auto rho = createTestFluidDensityModel(xml);

    ArrayType vars = {{473.15, 4.e+7}};
    const double rho_expected = 890.943136237744;
    ASSERT_NEAR(rho_expected, rho->getValue(vars), 1.e-10);

    const double drho_dT = rho->getdValue(vars, PropertyVariableType::T);
    const double drho_dp = rho->getdValue(vars, PropertyVariableType::p);

    const double perturbation = 1.e-4;

    // Test the differentiation: with respect to temperature:
    vars[static_cast<unsigned>(PropertyVariableType::T)] += perturbation;
    const double rho_T1 = rho->getValue(vars);
    ASSERT_NEAR((rho_T1 - rho_expected) / perturbation, drho_dT, 1.e-6);

    // Test the differentiation: with respect to pressure:
    vars[static_cast<unsigned>(PropertyVariableType::p)] += perturbation;
    const double rho_p1 = rho->getValue(vars);
    ASSERT_NEAR((rho_p1 - rho_T1) / perturbation, drho_dp, 1.e-6);
}
