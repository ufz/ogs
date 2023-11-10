/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialSpatialDistributionMap.h"
#include "Medium.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "Phase.h"

namespace MaterialPropertyLib
{
template <typename ContainerMedium,
          typename ContainerSolid,
          typename ContainerLiquid,
          typename ContainerGas>
void checkMaterialSpatialDistributionMap(
    MeshLib::Mesh const& mesh,
    MaterialPropertyLib::MaterialSpatialDistributionMap const& media_map,
    ContainerMedium const& required_properties_medium,
    ContainerSolid const& required_properties_solid_phase,
    ContainerLiquid const& required_properties_liquid_phase,
    ContainerGas const& required_properties_gas_phase)
{
    for (auto const element_id : mesh.getElements() | MeshLib::views::ids)
    {
        auto const& medium = *media_map.getMedium(element_id);
        if (!required_properties_medium.empty())
        {
            MaterialPropertyLib::checkRequiredProperties(
                medium, required_properties_medium);
        }
        if (!required_properties_liquid_phase.empty())
        {
            MaterialPropertyLib::checkRequiredProperties(
                medium.phase("AqueousLiquid"),
                required_properties_liquid_phase);
        }
        if (!required_properties_gas_phase.empty())
        {
            MaterialPropertyLib::checkRequiredProperties(
                medium.phase("Gas"), required_properties_gas_phase);
        }
        if (!required_properties_solid_phase.empty())
        {
            MaterialPropertyLib::checkRequiredProperties(
                medium.phase("Solid"), required_properties_solid_phase);
        }
    }
}

}  // namespace MaterialPropertyLib
