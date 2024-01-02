/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "BaseLib/Error.h"
#include "MeshLib/Location.h"
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
                       std::vector<T> const& values)
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
        name, item_type, number_of_components);
    if (!property)
    {
        OGS_FATAL("Error while creating PropertyVector '{:s}'.", name);
    }
    property->reserve(values.size());
    std::copy(values.cbegin(), values.cend(), std::back_inserter(*property));
}
}  // namespace MeshLib
