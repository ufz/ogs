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

#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/TwoPhaseModels/CreateTwoPhaseFlowMaterialProperties.h"

#include "ThermalTwoPhaseFlowWithPPMaterialProperties.h"
namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
namespace ThermalTwoPhaseFlowWithPP
{
std::unique_ptr<ThermalTwoPhaseFlowWithPPMaterialProperties>
createThermalTwoPhaseFlowWithPPMaterialProperties(
    BaseLib::ConfigTree const& config,
    MeshLib::PropertyVector<int> const& material_ids);

}  // end namespace
}  // end namespace
