/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   CreateLiquidFlowMaterialProperties.h
 *
 *   Created on December 14, 2016, 1:20 PM
 */


#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib
{
namespace LiquidFlow
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
    BaseLib::ConfigTree const& config, bool const has_material_ids,
    MeshLib::PropertyVector<int> const& material_ids);

}  // end of namespace
}  // end of namespace
