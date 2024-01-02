/**
 * \file
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

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/Fluid/ConstantFluidProperty.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/SpecificHeatCapacity/CreateSpecificFluidHeatCapacityModel.h"
#include "MaterialLib/PhysicalConstant.h"
#include "Tests/TestTools.h"

using namespace MaterialLib;
using namespace MaterialLib::Fluid;

using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

std::unique_ptr<FluidProperty> createSpecificFluidHeatCapacityModel(
    const char xml[])
{
    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree conf(std::move(ptree), "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("specific_heat_capacity");
    return createSpecificFluidHeatCapacityModel(sub_config);
}

TEST(MaterialFluid, checkConstantSpecificFluidHeatCapacityModel)
{
    const char xml[] =
        "<specific_heat_capacity>"
        "   <type>Constant</type>"
        "   <value> 900. </value> "
        "</specific_heat_capacity>";
    const auto cp = createSpecificFluidHeatCapacityModel(xml);

    ArrayType dummy;
    ASSERT_EQ(900., cp->getValue(dummy));
}
