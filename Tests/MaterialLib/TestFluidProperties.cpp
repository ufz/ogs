/**
 * \file
 *  \brief Test the creator of FluidProperties
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

#include <cmath>
#include <memory>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/Fluid/FluidProperties/CreateFluidProperties.h"
#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"
#include "Tests/TestTools.h"

using namespace MaterialLib;
using namespace MaterialLib::Fluid;

using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

std::unique_ptr<FluidProperties> createTestFluidProperties(const char xml[])
{
    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree conf(std::move(ptree), "", BaseLib::ConfigTree::onerror,
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
    const double dmu_dT_expected =
        -1.e-3 / 368 * std::exp(-(vars[0] - 293) / 368);
    ASSERT_NEAR(
        dmu_dT_expected,
        fluid_model->getdValue(FluidPropertyType::Viscosity, vars,
                               MaterialLib::Fluid::PropertyVariableType::T),
        1.e-10);

    vars[0] = 273.1;
    ASSERT_NEAR(1000.0 * (1 - 4.3e-4 * (vars[0] - 293.0)),
                fluid_model->getValue(FluidPropertyType::Density, vars),
                1.e-10);
    ASSERT_NEAR(-1000.0 * 4.3e-4,
                fluid_model->getdValue(FluidPropertyType::Density, vars,
                                       Fluid::PropertyVariableType::T),
                1.e-10);
}
