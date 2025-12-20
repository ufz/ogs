// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <vector>

#include "BaseLib/Error.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Properties.h"

namespace MeshLib
{
/// Creates a new \c PropertyVector in the given mesh and initializes it with
/// the given data. A \c PropertyVector with the same name must not exist.
/// \param mesh A \c Mesh the new \c ProperyVector will be created in.
/// \param name A string that contains the name of the new \c PropertyVector.
/// \param item_type One of the values \c MeshLib::MeshItemType::Cell or \c
/// \c MeshLib::MeshItemType::Node that shows the association of the property
/// values either to \c Element's / cells or \c Node's
/// \param number_of_components the number of components of a property
/// \param values A vector containing the values that are used for
/// initialization.
template <typename T>
void addPropertyToMesh(Mesh& mesh, std::string_view name,
                       MeshItemType item_type, std::size_t number_of_components,
                       std::span<T const> values)
{
    if (item_type == MeshItemType::Node)
    {
        if (mesh.getNumberOfNodes() != values.size() / number_of_components)
        {
            OGS_FATAL(
                "Error number of nodes ({:d}) does not match the number of "
                "tuples ({:d}).",
                mesh.getNumberOfNodes(), values.size() / number_of_components);
        }
    }
    if (item_type == MeshItemType::Cell)
    {
        if (mesh.getNumberOfElements() != values.size() / number_of_components)
        {
            OGS_FATAL(
                "Error number of elements ({:d}) does not match the number of "
                "tuples ({:d}).",
                mesh.getNumberOfElements(),
                values.size() / number_of_components);
        }
    }

    auto* const property = mesh.getProperties().createNewPropertyVector<T>(
        name, item_type, values.size() / number_of_components,
        number_of_components);
    if (!property)
    {
        OGS_FATAL("Error while creating PropertyVector '{:s}'.", name);
    }
    assert(property->size() == values.size());
    property->assign(values);
}
}  // namespace MeshLib
