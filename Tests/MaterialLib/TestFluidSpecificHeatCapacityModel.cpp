/**
 *  \copyright
 *   Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *   \file TestFluidSpecificHeatCapacityModel.cpp
 *
 */

#include <gtest/gtest.h>

#include <memory>

#include "Tests/TestTools.h"

#include "BaseLib/ConfigTree.h"

#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/ConstantFluidProperty.h"
#include "MaterialLib/Fluid/SpecificHeatCapacity/CreateSpecificFluidHeatCapacityModel.h"

using namespace MaterialLib;
using namespace MaterialLib::Fluid;

using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

std::unique_ptr<FluidProperty> createSpecificFluidHeatCapacityModel(
    const char xml[])
{
    auto const ptree = readXml(xml);
    BaseLib::ConfigTree conf(ptree, "", BaseLib::ConfigTree::onerror,
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
