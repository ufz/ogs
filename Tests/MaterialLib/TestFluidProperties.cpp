/**
 *  \brief Test the creator of FluidProperties
 *
 *  \copyright
 *   Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file  TestFluidProperties.cpp
 *
 */

#include <gtest/gtest.h>

#include <cmath>
#include <memory>

#include "Tests/TestTools.h"

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/Fluid/FluidProperties/CreateFluidProperties.h"
#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"

using namespace MaterialLib;
using namespace MaterialLib::Fluid;

using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

std::unique_ptr<FluidProperties> createTestFluidProperties(const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("fluid");
    return createFluidProperties(sub_config);
}

TEST(MaterialFluidProperties, checkPrimaryVariableDependentFluidProperties)
{
    const char xml[] =
        "<fluid>"
        "    <density>"
        "       <type>TemperatureDependent</type>"
        "       <temperature0> 293.0 </temperature0> "
        "       <beta> 4.3e-4 </beta> "
        "       <rho0>1000.</rho0>"
        "    </density>"
        "    <viscosity>"
        "       <type>TemperatureDependent</type>"
        "       <mu0>1.e-3 </mu0>"
        "       <tc>293.</tc>"
        "       <tv>368.</tv>"
        "    </viscosity>"
        "</fluid>";

    auto const fluid_model = createTestFluidProperties(xml);

    ArrayType vars;
    vars[0] = 350.0;
    const double mu_expected = 1.e-3 * std::exp(-(vars[0] - 293) / 368);
    ASSERT_NEAR(mu_expected,
                fluid_model->getValue(FluidPropertyType::Viscosity, vars),
                1.e-10);
    ASSERT_NEAR(
        -mu_expected,
        fluid_model->getdValue(FluidPropertyType::Viscosity, vars,
                               MaterialLib::Fluid::PropertyVariableType::T),
        1.e-10);

    vars[0] = 273.1;
    ASSERT_NEAR(1000.0 * (1 + 4.3e-4 * (vars[0] - 293.0)),
                fluid_model->getValue(FluidPropertyType::Density, vars),
                1.e-10);
    ASSERT_NEAR(1000.0 * 4.3e-4,
                fluid_model->getdValue(FluidPropertyType::Density, vars,
                                       Fluid::PropertyVariableType::T),
                1.e-10);
}

TEST(MaterialFluidProperties, checkFluidPropertiesWithDensityDependentModels_T)
{
    const char xml[] =
        "<fluid>"
        "    <density>"
        "        <type>TemperatureDependent</type>"
        "        <temperature0> 293.0 </temperature0> "
        "        <beta> 2.5003219164466073e-05 </beta> "
        "        <rho0> 998.</rho0>"
        "    </density>"
        "    <viscosity>"
        "       <type>WaterViscosityIAPWS</type>"
        "    </viscosity>"
        "</fluid>";

    auto const fluid_model = createTestFluidProperties(xml);

    ArrayType vars;
    vars[0] = 373.15;

    const double rho = fluid_model->getValue(FluidPropertyType::Density, vars);
    ASSERT_NEAR(1000.0, rho, 1.e-10);

    const double mu = fluid_model->getValue(FluidPropertyType::Viscosity, vars);
    ASSERT_NEAR(0.307883622e-3, mu, 1.e-10);

    const double drho_dT =
        fluid_model->getdValue(FluidPropertyType::Density, vars,
                               MaterialLib::Fluid::PropertyVariableType::T);
    const double dmu_dT =
        fluid_model->getdValue(FluidPropertyType::Viscosity, vars,
                               MaterialLib::Fluid::PropertyVariableType::T);

    const double perturbation = 1.e-6;
    vars[0] += perturbation;
    const double rho1 = fluid_model->getValue(FluidPropertyType::Density, vars);
    const double mu1 = fluid_model->getValue(FluidPropertyType::Viscosity, vars);

    ASSERT_NEAR((rho1 - rho) / perturbation, drho_dT, 1.e-7);
    ASSERT_NEAR((mu1 - mu) / perturbation, dmu_dT, 1.e-10);
}

TEST(MaterialFluidProperties, checkFluidPropertiesWithDensityDependentModels_dp)
{
    const char xml[] =
        "<fluid>"
        "   <density>"
        "       <type>LiquidDensity</type>"
        "       <temperature0> 273.15 </temperature0> "
        "       <p0> 1.e+5 </p0> "
        "       <bulk_modulus> 2.15e+9 </bulk_modulus> "
        "       <beta> 2.0e-4 </beta> "
        "       <rho0>999.8</rho0>"
        "    </density>"
        "    <viscosity>"
        "       <type>WaterViscosityIAPWS</type>"
        "    </viscosity>"
        "</fluid>";

    auto const fluid_model = createTestFluidProperties(xml);

    ArrayType vars = {{273.15 + 60.0, 1.e+6}};

    const double rho = fluid_model->getValue(FluidPropertyType::Density, vars);
    const double mu = fluid_model->getValue(FluidPropertyType::Viscosity, vars);
    const double drho_dp =
        fluid_model->getdValue(FluidPropertyType::Density, vars,
                               MaterialLib::Fluid::PropertyVariableType::p);
    const double dmu_dp =
        fluid_model->getdValue(FluidPropertyType::Viscosity, vars,
                               MaterialLib::Fluid::PropertyVariableType::p);

    const double perturbation = 1.e-6;
    vars[1] += perturbation;
    const double rho1 = fluid_model->getValue(FluidPropertyType::Density, vars);
    const double mu1 = fluid_model->getValue(FluidPropertyType::Viscosity, vars);

    // Only check d()/drho * drho/dp (FluidPropertiesWithDensityDependentModels)
    // The other functionalities of fluid property models are checked in other
    // tests in TestFluid*.cpp files.
    ASSERT_NEAR((rho1 - rho) / perturbation, drho_dp, 1.e-7);
    ASSERT_NEAR((mu1 - mu) / perturbation, dmu_dp, 1.e-10);
}
