// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once
#include <string>

#include "MeshLib/Mesh.h"
#include "MeshLib/PropertyVector.h"

namespace MeshLib
{
/// \returns a PropertyVector of the corresponding type, name on nodes, or
/// cells, or integration points if such exists, or creates a new one.
/// \attention For the integration points the result's size is zero.
/// \see MeshLib::addPropertyToMesh()
template <typename T>
PropertyVector<T>* getOrCreateMeshProperty(Mesh& mesh,
                                           std::string const& property_name,
                                           MeshItemType const item_type,
                                           int const number_of_components)
{
    if (property_name.empty())
    {
        OGS_FATAL(
            "Trying to get or to create a mesh property with empty name.");
    }

    auto numberOfMeshItems = [&mesh, &item_type]() -> std::size_t
    {
        switch (item_type)
        {
            case MeshItemType::Cell:
                return mesh.getNumberOfElements();
            case MeshItemType::Node:
                return mesh.getNumberOfNodes();
            case MeshItemType::IntegrationPoint:
                return 0;  // For the integration point data the size is
                           // variable
            default:
                OGS_FATAL(
                    "getOrCreateMeshProperty cannot handle other "
                    "types than Node, Cell, or IntegrationPoint.");
        }
        return 0;
    };

    if (mesh.getProperties().existsPropertyVector<T>(property_name))
    {
        auto result =
            mesh.getProperties().template getPropertyVector<T>(property_name);
        assert(result);
        if (item_type != MeshItemType::IntegrationPoint)
        {
            // Test the size if number of mesh items is known, which is not the
            // case for the integration point data.
            assert(result->size() ==
                   numberOfMeshItems() * number_of_components);
        }
        return result;
    }

    auto result = mesh.getProperties().template createNewPropertyVector<T>(
        property_name, item_type, numberOfMeshItems(), number_of_components);
    assert(result);
    return result;
}
}  // namespace MeshLib
