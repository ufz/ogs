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
        DBUG("Reusing existing property vector '{}' on mesh '{}'.",
             property_name, mesh.getName());

        auto* const result =
            mesh.getProperties().template getPropertyVector<T>(property_name);
        assert(result);
        if (auto const mit = result->getMeshItemType(); mit != item_type)
        {
            OGS_FATAL(
                "Found mesh item type '{}' for mesh property '{}' on mesh "
                "'{}', but expected a '{}' property.",
                toString(mit), property_name, mesh.getName(),
                toString(item_type));
        }
        if (auto const ncomp = result->getNumberOfGlobalComponents();
            ncomp != number_of_components)
        {
            OGS_FATAL(
                "Found {} components for mesh property '{}' on mesh '{}', but "
                "expected {} components.",
                ncomp, property_name, mesh.getName(), number_of_components);
        }
        if (item_type != MeshItemType::IntegrationPoint)
        {
            // Test the size if number of mesh items is known, which is not the
            // case for the integration point data.
            auto const size = result->size();
            auto const size_expected =
                numberOfMeshItems() * number_of_components;
            if (size != size_expected)
            {
                OGS_FATAL(
                    "Actual and expected size of property '{}' on mesh '{}' do "
                    "not match: {} != {}.",
                    property_name, mesh.getName(), size, size_expected);
            }
        }
        return result;
    }

    auto* const result =
        mesh.getProperties().template createNewPropertyVector<T>(
            property_name, item_type, numberOfMeshItems(),
            number_of_components);
    assert(result);
    return result;
}
}  // namespace MeshLib
