/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file
 *
 *   Created on December 14, 2016, 1:20 PM
 */

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ParameterLib
{
struct ParameterBase;
}

namespace ProcessLib::LiquidFlow
{
class LiquidFlowMaterialProperties;

/**
 *  Parse the XML input for fluid properties of a single phase and create an
 *  instance of class LiquidFlowMaterialProperties.
 *
 *  The XML syntax example is given in an unit test in
 *  Tests/Process/LiquidFlow/TestLiquidFlowMaterialProperties.cpp
 */
std::unique_ptr<LiquidFlowMaterialProperties>
createLiquidFlowMaterialProperties(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const& parameters,
    MeshLib::PropertyVector<int> const* const material_ids);

}  // namespace ProcessLib::LiquidFlow
