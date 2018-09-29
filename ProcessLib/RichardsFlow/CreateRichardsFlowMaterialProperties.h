/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "RichardsFlowMaterialProperties.h"
namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
namespace RichardsFlow
{
std::unique_ptr<RichardsFlowMaterialProperties>
createRichardsFlowMaterialProperties(
    BaseLib::ConfigTree const& config,
    MeshLib::PropertyVector<int> const* material_ids,
    std::vector<std::unique_ptr<ParameterBase>> const& parameters);

}  // end namespace
}  // end namespace
