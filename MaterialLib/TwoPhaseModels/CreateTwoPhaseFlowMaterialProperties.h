/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <tuple>

#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"

#include "TwoPhaseFlowWithPPMaterialProperties.h"
namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace TwoPhaseFlowWithPP
{
std::tuple<std::unique_ptr<TwoPhaseFlowWithPPMaterialProperties>,
           BaseLib::ConfigTree>
createTwoPhaseFlowMaterialProperties(
    BaseLib::ConfigTree const& config,
    MeshLib::PropertyVector<int> const& material_ids,
    std::vector<std::unique_ptr<ProcessLib::ParameterBase>> const& parameters);

}  // end namespace
}  // end namespace
