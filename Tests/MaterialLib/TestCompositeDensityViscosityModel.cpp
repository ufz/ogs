/**
 *  \brief Test Composite density viscosity models
 *
 *  \copyright
 *   Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file  TestCompositeDensityViscosityModel.cpp
 *
 */

#include <gtest/gtest.h>

#include <memory>
#include <cmath>

#include "TestTools.h"

#include "MaterialLib/Fluid/Density/createFluidDensityModel.h"
#include "MaterialLib/Fluid/Viscosity/createViscosityModel.h"
#include "MaterialLib/Fluid/CompositeFluidProperty.h"
#include "MaterialLib/Fluid/CompositeFluidProperty/CompositeDensityViscosityModel.h"

using namespace MaterialLib;
using namespace MaterialLib::Fluid;
using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

template <typename F>
std::unique_ptr<FluidProperty> createTestModel(const char xml[], F func,
                                               const std::string& key)
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree(key);
    return func(sub_config);
}

TEST(MaterialFluidModel, checkCompositeDensityViscosityModel)
{
    const char xml_d[] =
        "<density>"
        "   <type>TemperatureDependent</type>"
        "   <temperature0> 293.0 </temperature0> "
        "   <beta> 4.3e-4 </beta> "
        "   <rho0>1000.</rho0>"
        "</density>";

    auto rho = createTestModel(xml_d, createFluidDensityModel, "density");

    const char xml_v[] =
        "<viscosity>"
        "  <type>TemperatureDependent</type>"
        "  <mu0>1.e-3 </mu0>"
        "   <tc>293.</tc>"
        "   <tv>368.</tv>"
        "</viscosity>";
    auto mu = createTestModel(xml_v, createViscosityModel, "viscosity");

    std::unique_ptr<CompositeFluidProperty> composite_fluid_model =
        std::unique_ptr<CompositeFluidProperty>(
            new CompositeDensityViscosityModel(std::move(rho), std::move(mu)));

    ArrayType vars;
    vars[0] = 350.0;
    const double mu_expected = 1.e-3 * std::exp(-(vars[0] - 293) / 368);
    ASSERT_NEAR(mu_expected,
                composite_fluid_model->getValue(PropertyType::Vicosity, vars),
                1.e-10);
    ASSERT_NEAR(-mu_expected,
                composite_fluid_model->getdValue(
                    PropertyType::Vicosity, vars,
                    MaterialLib::Fluid::PropertyVariableType::T),
                1.e-10);

    vars[0] = 273.1;
    ASSERT_NEAR(1000.0 * (1 + 4.3e-4 * (vars[0] - 293.0)),
                composite_fluid_model->getValue(PropertyType::Density, vars),
                1.e-10);
    ASSERT_NEAR(1000.0 * 4.3e-4, composite_fluid_model->getdValue(
                                     PropertyType::Density, vars,
                                     Fluid::PropertyVariableType::T),
                1.e-10);
}
