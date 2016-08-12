/*!
   \file  TestFluidProperties.cpp
   \brief Test classes for fluid properties.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#include <gtest/gtest.h>

#include <memory>
#include <algorithm>  // std::max

#include "TestTools.h"

#include "MaterialLib/Fluid/Density/createFluidDensityModel.h"
#include "MaterialLib/Fluid/Viscosity/createViscosityModel.h"
#include "MaterialLib/PhysicalConstant.h"

namespace
{
using namespace MaterialLib;
using namespace MaterialLib::Fluid;

//----------------------------------------------------------------------------
// Test density models.
FluidProperty* createTestFluidDensityModel(const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("fluid_density");
    return MaterialLib::Fluid::createFluidDensityModel(&sub_config);
}

TEST(Material, checkConstantDensity)
{
    const char xml[] =
        "<fluid_density>"
        "   <type>constant</type>"
        "   <value> 998.0 </value> "
        "</fluid_density>";
    std::unique_ptr<FluidProperty> rho(createTestFluidDensityModel(xml));

    ASSERT_EQ(998.0, rho->getValue(nullptr));
    ASSERT_EQ(0.0, rho->getdValue(nullptr,
                                  MaterialLib::Fluid::PropertyVariable::T));
}

TEST(Material, checkIdealGasLaw)
{
    const char xml[] =
        "<fluid_density>"
        "   <type>ideal_gas_law</type>"
        "   <molar_mass> 28.96 </molar_mass> "
        "</fluid_density>";
    std::unique_ptr<FluidProperty> rho(createTestFluidDensityModel(xml));

    const double molar_air = 28.96;
    const double T = 290.;
    const double p = 1.e+5;
    const double R = PhysicalConstant::IdealGasConstant;
    const double expected_air_dens = molar_air * p / (R * T);
    double vars[] = {290, 0, 1.e+5};
    ASSERT_NEAR(expected_air_dens, rho->getValue(vars), 1.e-10);

    const double expected_d_air_dens_dT = -molar_air * p / (R * T * T);
    ASSERT_NEAR(expected_d_air_dens_dT, rho->getdValue(vars,
                                        Fluid::PropertyVariable::T), 1.e-10);

    const double expected_d_air_dens_dp = molar_air / (R * T);
    ASSERT_NEAR(expected_d_air_dens_dp,
                rho->getdValue(vars, Fluid::PropertyVariable::pg), 1.e-10);
}

TEST(Material, checkLinearTemperatureDependentDensity)
{
    const char xml[] =
        "<fluid_density>"
        "   <type>temperature_dependent</type>"
        "   <temperature0> 293.0 </temperature0> "
        "   <beta> 4.3e-4 </beta> "
        "   <rho0>1000.</rho0>"
        "</fluid_density>";

    std::unique_ptr<FluidProperty> rho(createTestFluidDensityModel(xml));

    double vars[] = {273.1 + 60.0};
    ASSERT_NEAR(1000.0 * (1 + 4.3e-4 * (vars[0] - 293.0)), rho->getValue(vars),
                1.e-10);
    ASSERT_NEAR(1000.0 * 4.3e-4, rho->getdValue(vars,
                                         Fluid::PropertyVariable::T), 1.e-10);
}

TEST(Material, checkLiquidDensity)
{
    const char xml[] =
        "<fluid_density>"
        "   <type>liquid_density</type>"
        "   <temperature0> 273.15 </temperature0> "
        "   <p0> 1.e+5 </p0> "
        "   <bulk_modulus> 2.15e+9 </bulk_modulus> "
        "   <beta> 2.0e-4 </beta> "
        "   <rho0>999.8</rho0>"
        "</fluid_density>";
    std::unique_ptr<FluidProperty> rho(createTestFluidDensityModel(xml));

    const double vars[] = {273.15 + 60.0, 1.e+6};
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
                rho->getdValue(vars, Fluid::PropertyVariable::T), 1.e-10);

    // Test the derivative with respect to pressure.
    const double fac_p = 1. - (p - p0) / K;
    ASSERT_NEAR(rho0 / (1. + beta * (T - T0)) / (fac_p * fac_p * K),
                rho->getdValue(vars, Fluid::PropertyVariable::pl), 1.e-10);    
}

//----------------------------------------------------------------------------
// Test viscosity models.
FluidProperty* createTestViscosityModel(const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("viscosity");
    return MaterialLib::Fluid::createViscosityModel(&sub_config);
}

TEST(Material, checkConstantViscosity)
{
    const char xml[] =
        "<viscosity>"
        "   <type>constant</type>"
        "   <value> 1.e-4 </value> "
        "</viscosity>";
    std::unique_ptr<FluidProperty> mu(createTestViscosityModel(xml));

    ASSERT_EQ(1.e-4, mu->getValue(nullptr));
    ASSERT_EQ(0.0, mu->getdValue(nullptr,
                                  MaterialLib::Fluid::PropertyVariable::T));
}

TEST(Material, checkTemperatureDependentViscosity)
{
    const char xml[] =
        "<viscosity>"
        "  <type>temperature_dependent</type>"
        "  <mu0>1.e-3 </mu0>"
        "   <tc>293.</tc>"
        "   <tv>368.</tv>"
        "</viscosity>";
    std::unique_ptr<FluidProperty> mu(createTestViscosityModel(xml));

    const double vars[] = {350.0};
    const double mu_expected = 1.e-3 * std::exp(-(vars[0] - 293) / 368);
    // Test the density.
    ASSERT_NEAR(mu_expected, mu->getValue(vars), 1.e-10);
    // Test the derivative with respect to temperature.
    ASSERT_NEAR(-mu_expected, mu->getdValue(vars,
                              MaterialLib::Fluid::PropertyVariable::T), 1.e-10);
}

TEST(Material, checkLinearPressureDependentViscosity)
{
    const char xml[] =
        "<viscosity>"
        "  <type>linear_pressure</type>"
        "  <mu0>1.e-3 </mu0>"
        "   <p0>1.e+5</p0>"
        "   <gamma>1.e-6</gamma>"
        "</viscosity>";
    std::unique_ptr<FluidProperty> mu(createTestViscosityModel(xml));

    const double vars[] = {293, 2.e+6};
    // Test the density.
    ASSERT_NEAR(1.e-3 * (1. + 1.e-6 * (vars[1] - 1.e+5)),
                                                   mu->getValue(vars), 1.e-10);
    // Test the derivative with respect to pressure.
    ASSERT_NEAR(1.e-9, mu->getdValue(vars,
                             MaterialLib::Fluid::PropertyVariable::pl), 1.e-10);
}

TEST(Material, checkVogelViscosity)
{
    const char xml_w[] =
        "<viscosity>"
        "  <type>vogels</type>"
        "  <liquid_type>water </liquid_type>"
        "</viscosity>";
    std::unique_ptr<FluidProperty> mu_w(createTestViscosityModel(xml_w));
    double vars[] = {303.0};
    const auto var_type = MaterialLib::Fluid::PropertyVariable::T;
    ASSERT_NEAR(0.802657e-3, mu_w->getValue(vars), 1.e-5);
    ASSERT_NEAR(-1.87823e-5, mu_w->getdValue(vars, var_type), 1.e-5);

    const char xml_co2[] =
        "<viscosity>"
        "  <type>vogels</type>"
        "  <liquid_type> CO2 </liquid_type>"
        "</viscosity>";
    std::unique_ptr<FluidProperty> mu_co2(createTestViscosityModel(xml_co2));
    vars[0] = 255.04; 
    ASSERT_NEAR(0.137956e-3, mu_co2->getValue(vars), 1.e-5);
    ASSERT_NEAR(-2.35664e-6, mu_co2->getdValue(vars, var_type), 1.e-5);

    const char xml_ch4[] =
        "<viscosity>"
        "  <type>vogels</type>"
        "  <liquid_type> CH4 </liquid_type>"
        "</viscosity>";
    std::unique_ptr<FluidProperty> mu_ch4(createTestViscosityModel(xml_ch4));
    vars[0] = 172.0; 
    ASSERT_NEAR(0.352072e-4, mu_ch4->getValue(vars), 1.e-5);
    ASSERT_NEAR(-2.35664e-6, mu_ch4->getdValue(vars, var_type), 1.e-5);
}

}
