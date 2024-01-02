/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 6, 2023, 4:21 PM
 */

#include "ZeroMeshFieldDataByMaterialIDs.h"

#include <algorithm>
#include <functional>

#include "IntegrationPointDataTools.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Utils/IntegrationPointWriter.h"

namespace MeshToolsLib
{
void zeroMeshFieldDataByMaterialIDs(
    MeshLib::Mesh& mesh, std::vector<int> const& selected_material_ids)
{
    auto const materialIds = materialIDs(mesh);
    if (!materialIds)
    {
        OGS_FATAL(
            "Mesh contains no int-property vector named 'MaterialIDs' needed "
            "by the 'zero_mesh_field_data_by_materialIDs'");
    }

    std::vector<std::size_t> element_ids_for_selected_materials;
    for (std::size_t i = 0; i < materialIds->size(); ++i)
    {
        if (std::find(selected_material_ids.begin(),
                      selected_material_ids.end(),
                      (*materialIds)[i]) != selected_material_ids.end())
        {
            element_ids_for_selected_materials.push_back(i);
        }
    }

    MeshLib::Properties& properties = mesh.getProperties();

    std::vector<std::size_t> element_ip_data_offsets;

    for (auto const& [name, property] : properties)
    {
        if (auto const item_type = property->getMeshItemType();
            item_type != MeshLib::MeshItemType::IntegrationPoint)
        {
            continue;
        }

        // For special field data such as OGS_VERSION,
        // IntegrationPointMetaData,
        // etc., which are not "real" integration points:
        if (!property->getPropertyName().ends_with("_ip"))
        {
            continue;
        }

        if (properties.template hasPropertyVector<double>(
                name, MeshLib::MeshItemType::IntegrationPoint))
        {
            auto& pv = *properties.template getPropertyVector<double>(name);
            const int n_components = pv.getNumberOfGlobalComponents();

            if (element_ip_data_offsets.empty())
            {
                // The returned values has already been multiplied with
                // pv.getNumberOfGlobalComponents()
                element_ip_data_offsets =
                    MeshToolsLib::getIntegrationPointDataOffsetsOfMeshElements(
                        mesh.getElements(), pv, properties);

                // element_ip_data_offsets / pv.getNumberOfGlobalComponents()
                std::transform(element_ip_data_offsets.begin(),
                               element_ip_data_offsets.end(),
                               element_ip_data_offsets.begin(),
                               [n = n_components](double const v)
                               { return v / n; });
            }

            for (auto const element_id : element_ids_for_selected_materials)
            {
                std::fill(
                    pv.begin() +
                        n_components * element_ip_data_offsets[element_id],
                    pv.begin() +
                        n_components * element_ip_data_offsets[element_id + 1],
                    0.0);
            }
        }
    }
}
}  // namespace MeshToolsLib
