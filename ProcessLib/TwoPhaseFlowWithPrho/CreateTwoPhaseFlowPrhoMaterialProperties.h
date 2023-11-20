/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "ProcessLib/TwoPhaseFlowWithPrho/TwoPhaseFlowWithPrhoMaterialProperties.h"
namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
namespace TwoPhaseFlowWithPrho
{
std::unique_ptr<TwoPhaseFlowWithPrhoMaterialProperties>
createTwoPhaseFlowPrhoMaterialProperties(
    BaseLib::ConfigTree const& config,
    MeshLib::PropertyVector<int> const* material_ids);

}  // namespace TwoPhaseFlowWithPrho
}  // namespace ProcessLib
