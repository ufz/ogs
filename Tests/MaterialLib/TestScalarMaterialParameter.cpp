/*!
   \file  TestScalarMaterialParameter.cpp
   \brief Test classes for scalar material parameters.

   \author Wenqing Wang
   \date Jan 2015

   \copyright
    Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
               Distributed under a Modified BSD License.
               See accompanying file LICENSE.txt or
               http://www.opengeosys.org/project/license
*/
#include <gtest/gtest.h>

#include <algorithm>  // std::max

#include "TestTools.h"

#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/ScalarParameter.h"
#include "MaterialLib/ConstantScalarModel.h"

#include "MaterialLib/Solid/Density/LinearSolidDensityModel.h"
#include "MaterialLib/Fluid/Density/IdealGasLaw.h"
#include "MaterialLib/Fluid/Density/LinearTemperatureDependentDensity.h"
#include "MaterialLib/Fluid/Density/LiquidDensity.h"

#include "MaterialLib/Fluid/Viscosity/LinearPressureDependentViscosity.h"
#include "MaterialLib/Fluid/Viscosity/TemperatureDependentViscosity.h"
#include "MaterialLib/Fluid/Viscosity/VogelsLiquidDynamicViscosity.h"

namespace
{
// For local assembly
template <typename T_MAT>
class TestScalar
{
public:
    TestScalar(const T_MAT& mat) : _mat(const_cast<T_MAT*>(&mat)) {}
    // Test function
    template <typename... Args>
    double getMatParameterTest(Args... args) const
    {
        return _mat->getValue(args...);
    }

private:
    T_MAT* _mat;
};

using namespace MaterialLib;
using namespace MaterialLib::Fluid;

TEST(Material, checkDensity)
{
    //-- Solid --------------------------------------------------------------
    // Constant
    const double rho = 2080.;
    MaterialLib::ScalarParameter<FluidDensityType, ConstantScalarModel>
        s_density(rho);
    ASSERT_NEAR(rho, s_density.getValue(), 1.e-10);
    ASSERT_EQ(FluidDensityType::CONSTANT, s_density.getType());

    // Linear
    MaterialLib::ScalarParameter<FluidDensityType, LinearSolidDensityModel>
        lin_density(20., 1000., 100., 900.);
    const double lin_den_expected =
        (900. - 1000.) / (100. - 20.) * (50. - 20.) + 1000.;
    ASSERT_NEAR(lin_den_expected, lin_density.getValue(50.), 1.e-10);

    //-- Fluid --------------------------------------------------------------
    const double molar_air = 28.96;
    MaterialLib::ScalarParameter<FluidDensityType, IdealGasLaw> air_density(
        molar_air);
    const double T = 290.;
    const double p = 1.e+5;
    const double R = PhysicalConstant::IdealGasConstant;
    const double expected_air_dens = molar_air * p / (R * T);
    ASSERT_NEAR(expected_air_dens, air_density.getValue(T, p), 1.e-10);

    const double expected_d_air_dens_dT = -molar_air * p / (R * T * T);
    ASSERT_NEAR(expected_d_air_dens_dT, air_density.getdValue(T, p, 0), 1.e-10);

    const double expected_d_air_dens_dp = molar_air / (R * T);
    ASSERT_NEAR(expected_d_air_dens_dp, air_density.getdValue(T, p, 1), 1.e-10);

    TestScalar<MaterialLib::ScalarParameter<FluidDensityType, IdealGasLaw>>
        test0(air_density);
    ASSERT_NEAR(expected_air_dens, test0.getMatParameterTest(T, p), 1.e-10);

    // Check polymorphy
    MaterialLib::ParameterBase<FluidDensityType>* den_ptr_base = &air_density;
    if (den_ptr_base->getType() == FluidDensityType::IDEAL_GAS)
    {
        MaterialLib::ScalarParameter<FluidDensityType, IdealGasLaw>* den_ptr =
            static_cast<
                MaterialLib::ScalarParameter<FluidDensityType, IdealGasLaw>*>(
                den_ptr_base);
        ASSERT_NEAR(molar_air * 2.e+5 / (R * 253.),
                    den_ptr->getValue(253., 2.e+5), 1.e-10);
    }
}

TEST(Material, checkLinearTemperatureDependentDensity)
{
    const char xml[] =
        "<fluid_density>"
        "   <type>linear_temperature</type>"
        "   <temperature0> 293.0 </temperature0> "
        "   <beta> 4.3e-4 </beta> "
        "   <rho0>1000.</rho0>"
        "</fluid_density>";
    auto const ptree = readXml(xml);

    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);

    auto const& sub_config = conf.getConfigSubtree("fluid_density");

    auto const type = sub_config.getConfigParameter<std::string>("type");
    const std::string exepected_type = "linear_temperature";
    ASSERT_STREQ(exepected_type.c_str(), type.c_str());

    MaterialLib::ScalarParameter<FluidDensityType,
                                 LinearTemperatureDependentDensity>
        density(&sub_config);
    const double T = 273.1 + 60.0;
    ASSERT_NEAR(1000.0 * (1 + 4.3e-4 * (T - 293.0)), density.getValue(T),
                1.e-10);
    ASSERT_NEAR(1000.0 * 4.3e-4, density.getdValue(T), 1.e-10);
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
    auto const ptree = readXml(xml);

    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);

    auto const& sub_config = conf.getConfigSubtree("fluid_density");

    auto const type = sub_config.getConfigParameter<std::string>("type");
    const std::string exepected_type = "liquid_density";
    ASSERT_STREQ(exepected_type.c_str(), type.c_str());

    MaterialLib::ScalarParameter<FluidDensityType, LiquidDensity> density(
        &sub_config);
    const double T = 273.15 + 60.0;
    const double p = 1.e+6;
    const double T0 = 273.15;
    const double p0 = 1.e+5;
    const double rho0 = 999.8;
    const double K = 2.15e+9;
    const double beta = 2.e-4;

    const double fac_T = 1. + beta * (T - T0);
    ASSERT_NEAR(rho0 / fac_T / (1. - (p - p0) / K), density.getValue(T, p),
                1.e-10);

    // Test the derivative with respect to temperature.
    ASSERT_NEAR(-beta * rho0 / (fac_T * fac_T) / (1. - (p - p0) / K),
                density.getdValue(T, p, 0), 1.e-10);

    // Test the derivative with respect to pressure.
    const double fac_p = 1. - (p - p0) / K;
    ASSERT_NEAR(rho0 / (1. + beta * (T - T0)) / (fac_p * fac_p * K),
                density.getdValue(T, p, 1), 1.e-10);
}

TEST(Material, checkViscosity)
{
    // Constant
    const double mu = 1.e-3;
    MaterialLib::ScalarParameter<ViscosityType, ConstantScalarModel> mu_const(
        mu);
    ASSERT_NEAR(mu, mu_const.getValue(), 1.e-10);
    ASSERT_EQ(ViscosityType::CONSTANT, mu_const.getType());

    const char xml[] =
        "<viscosity>"
        "  <type>temperature_dependent</type>"
        "  <mu0>1.e-3 </mu0>"
        "   <tc>293.</tc>"
        "   <tv>368.</tv>"
        "</viscosity>";
    auto const ptree = readXml(xml);

    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);

    auto const& sub_config = conf.getConfigSubtree("viscosity");
    auto const type = sub_config.getConfigParameter<std::string>("type");
    const std::string exepected_type = "temperature_dependent";
    ASSERT_STREQ(exepected_type.c_str(), type.c_str());

    MaterialLib::ScalarParameter<ViscosityType, TemperatureDependentViscosity>
        vis(&sub_config);

    const double T = 350.0;
    const double mu_expected = mu * std::exp(-(T - 293) / 368);
    // Test the density.
    ASSERT_NEAR(mu_expected, vis.getValue(T), 1.e-10);
    // Test the derivative with respect to temperature.
    ASSERT_NEAR(-mu_expected, vis.getdValue(T), 1.e-10);
}

TEST(Material, checkLinearPressureDependentViscosity)
{
    // Constant
    const double mu = 1.e-3;
    MaterialLib::ScalarParameter<ViscosityType, ConstantScalarModel> mu_const(
        mu);
    ASSERT_NEAR(mu, mu_const.getValue(), 1.e-10);
    ASSERT_EQ(ViscosityType::CONSTANT, mu_const.getType());

    const char xml[] =
        "<viscosity>"
        "  <type>linear_pressure</type>"
        "  <mu0>1.e-3 </mu0>"
        "   <p0>1.e+5</p0>"
        "   <gamma>1.e-6</gamma>"
        "</viscosity>";
    auto const ptree = readXml(xml);

    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);

    auto const& sub_config = conf.getConfigSubtree("viscosity");
    auto const type = sub_config.getConfigParameter<std::string>("type");
    const std::string exepected_type = "linear_pressure";
    ASSERT_STREQ(exepected_type.c_str(), type.c_str());

    MaterialLib::ScalarParameter<ViscosityType,
                                 LinearPressureDependentViscosity>
        vis(&sub_config);

    const double p = 2.e+6;
    // Test the density.
    ASSERT_NEAR(1.e-3 * (1. + 1.e-6 * (p - 1.e+5)), vis.getValue(p), 1.e-10);
    // Test the derivative with respect to pressure.
    ASSERT_NEAR(1.e-9, vis.getdValue(p), 1.e-10);
}

TEST(Material, checkVogelViscosity)
{
    VogelsLiquidDynamicViscosity mu_water(0);
    ASSERT_NEAR(0.802657e-3, mu_water.getValue(303.0), 1.e-5);
    ASSERT_NEAR(-1.87823e-5, mu_water.getdValue(303.0), 1.e-5);

    VogelsLiquidDynamicViscosity mu_CO2(1);
    ASSERT_NEAR(0.137956e-3, mu_CO2.getValue(255.04), 1.e-5);
    ASSERT_NEAR(-2.35664e-6, mu_CO2.getdValue(255.04), 1.e-5);

    VogelsLiquidDynamicViscosity mu_CH4(2);
    ASSERT_NEAR(0.352072e-4, mu_CH4.getValue(172), 1.e-5);
    ASSERT_NEAR(-2.35664e-6, mu_CH4.getdValue(172), 1.e-5);
}
}
