/**
 *  \brief Test viscosity models
 *
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file  TestFluidViscosity.cpp
 *
 */

#include <gtest/gtest.h>

#include <memory>
#include <cmath>

#include "TestTools.h"

#include "MaterialLib/Fluid/Viscosity/createViscosityModel.h"
#include "MaterialLib/PhysicalConstant.h"

using namespace MaterialLib;
using namespace MaterialLib::Fluid;
using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

std::unique_ptr<FluidProperty> createTestViscosityModel(const char xml[])
{
    auto const ptree = readXml(xml);
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
    ASSERT_EQ(0.0,
              mu->getdValue(dummy, MaterialLib::Fluid::PropertyVariableType::T));
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
    ASSERT_NEAR(-mu_expected,
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
    ASSERT_NEAR(1.e-9,
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
