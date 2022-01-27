/**
 *  \copyright
 *   Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *   \file
 *
 */

#include <gtest/gtest.h>

#include <memory>

#include "BaseLib/ConfigTree.h"
#include "MaterialLib/Fluid/ConstantFluidProperty.h"
#include "MaterialLib/Fluid/Density/LiquidDensity.h"
#include "MaterialLib/Fluid/FluidProperty.h"
#include "MaterialLib/Fluid/ThermalConductivity/CreateFluidThermalConductivityModel.h"
#include "MaterialLib/PhysicalConstant.h"
#include "Tests/TestTools.h"

using namespace MaterialLib;
using namespace MaterialLib::Fluid;

using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

std::unique_ptr<FluidProperty> createFluidThermalConductivityModel(
    const char xml[])
{
    auto ptree = Tests::readXml(xml);
    BaseLib::ConfigTree conf(std::move(ptree), "", BaseLib::ConfigTree::onerror,
                             BaseLib::ConfigTree::onwarning);
    auto const& sub_config = conf.getConfigSubtree("thermal_conductivity");
    return createFluidThermalConductivityModel(sub_config);
}

TEST(Material, checkConstantFluidThermalConductivity)
{
    const char xml[] =
        "<thermal_conductivity>"
        "   <type>Constant</type>"
        "   <value> .45 </value> "
        "</thermal_conductivity>";
    const auto lambda = createFluidThermalConductivityModel(xml);

    ArrayType dummy;
    ASSERT_EQ(.45, lambda->getValue(dummy));
}
