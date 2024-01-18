/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string_view>

#include "BaseLib/Error.h"
#include "MeshLib/Properties.h"
#include "MeshLib/Utils/IntegrationPointWriter.h"

namespace ProcessLib
{
/// Removes the suffix '_ip' from the passed field name.
inline std::string_view removeIPFieldDataNameSuffix(std::string_view const name)
{
    if (!name.ends_with("_ip"))
    {
        OGS_FATAL(
            "The name of integration point data must end with '_ip'. '{}' "
            "does not.",
            name);
    }

    return {name.data(), name.size() - 3};
}

template <typename LocalAssemblersVector>
void setIPDataInitialConditions(
    std::vector<std::unique_ptr<MeshLib::IntegrationPointWriter>> const&
        _integration_point_writer,
    MeshLib::Properties const& mesh_properties,
    LocalAssemblersVector& local_assemblers)
{
    for (auto const& ip_writer : _integration_point_writer)
    {
        // Find the mesh property with integration point writer's name.
        auto const& name = ip_writer->name();
        if (!mesh_properties.existsPropertyVector<double>(name))
        {
            continue;
        }
        auto const& mesh_property =
            *mesh_properties.template getPropertyVector<double>(name);

        // The mesh property must be defined on integration points.
        if (mesh_property.getMeshItemType() !=
            MeshLib::MeshItemType::IntegrationPoint)
        {
            continue;
        }

        auto const ip_meta_data =
            getIntegrationPointMetaData(mesh_properties, name);

        // Check the number of components.
        if (ip_meta_data.n_components !=
            mesh_property.getNumberOfGlobalComponents())
        {
            OGS_FATAL(
                "Different number of components in meta data ({:d}) than in "
                "the integration point field data for '{:s}': {:d}.",
                ip_meta_data.n_components, name,
                mesh_property.getNumberOfGlobalComponents());
        }

        INFO("Setting initial integration point data for '{}'", name);

        auto const name_transformed = removeIPFieldDataNameSuffix(name);

        // Now we have a properly named vtk's field data array and the
        // corresponding meta data.
        std::size_t position = 0;
        for (auto& local_asm : local_assemblers)
        {
            std::size_t const integration_points_read =
                local_asm->setIPDataInitialConditions(
                    name_transformed, &mesh_property[position],
                    ip_meta_data.integration_order);
            // The number of read integration points could be zero in case there
            // are e.g. multiple materials with different sets of internal state
            // variables.

            position += integration_points_read * ip_meta_data.n_components;
        }
    }
}
}  // namespace ProcessLib
