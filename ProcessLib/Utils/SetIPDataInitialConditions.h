/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "BaseLib/Error.h"
#include "MeshLib/Properties.h"
#include "MeshLib/Utils/IntegrationPointWriter.h"

namespace ProcessLib
{
template <typename LocalAssemblersVector>
void setIPDataInitialConditions(
    std::vector<std::unique_ptr<IntegrationPointWriter>> const&
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

        // Now we have a properly named vtk's field data array and the
        // corresponding meta data.
        std::size_t position = 0;
        for (auto& local_asm : local_assemblers)
        {
            std::size_t const integration_points_read =
                local_asm->setIPDataInitialConditions(
                    name, &mesh_property[position],
                    ip_meta_data.integration_order);
            if (integration_points_read == 0)
            {
                OGS_FATAL(
                    "No integration points read in the integration point "
                    "initial conditions set function.");
            }
            position += integration_points_read * ip_meta_data.n_components;
        }
    }
}
}  // namespace ProcessLib
