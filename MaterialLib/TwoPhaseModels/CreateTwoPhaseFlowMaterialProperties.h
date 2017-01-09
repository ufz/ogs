/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef OGS_CREATETWOPHASEFLOWMATERIALPROPERTIES_H
#define OGS_CREATETWOPHASEFLOWMATERIALPROPERTIES_H

#include <memory>
#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/TwoPhaseModels/TwoPhaseFlowWithPPMaterialProperties.h"
namespace BaseLib
{
class ConfigTree;
}

namespace MaterialLib
{
namespace TwoPhaseFlowWithPP
{
std::unique_ptr<TwoPhaseFlowWithPPMaterialProperties>
CreateTwoPhaseFlowMaterialProperties(
    BaseLib::ConfigTree const& config,
    bool const has_material_ids,
    MeshLib::PropertyVector<int> const& material_ids);

}  // end namespace
}  // end namespace

#endif /* CREATETWOPHASEFLOWMATERIALPROPERTIES_H */
