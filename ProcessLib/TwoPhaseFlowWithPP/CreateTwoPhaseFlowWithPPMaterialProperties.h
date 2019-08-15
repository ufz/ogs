/**
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "TwoPhaseFlowWithPPMaterialProperties.h"
namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
namespace TwoPhaseFlowWithPP
{
std::unique_ptr<TwoPhaseFlowWithPPMaterialProperties>
createTwoPhaseFlowWithPPMaterialProperties(
    BaseLib::ConfigTree const& config,
    MeshLib::PropertyVector<int> const* const material_ids,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);

}  // namespace TwoPhaseFlowWithPP
}  // namespace ProcessLib
