/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-28
 * \brief  Definition of the MeshInformation class
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>
#include <limits>
#include <map>
#include <optional>
#include <string>

#include "GeoLib/AABB.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Properties.h"

namespace MeshLib
{
/**
 * \brief A set of tools for extracting information from a mesh
 */
class MeshInformation
{
public:
    /// Returns the smallest and largest value of a scalar array with the given
    /// name.
    template <typename T>
    static std::optional<std::pair<T, T>> const getValueBounds(
        PropertyVector<T> const& property)
    {
        if (property.empty())
        {
            INFO("Mesh property vector '{:s}' is empty.",
                 property.getPropertyName());
            return std::nullopt;
        }

        auto const [min, max] =
            std::minmax_element(begin(property), end(property));
        return {{*min, *max}};
    }

    /// Returns the bounding box of the mesh.
    static GeoLib::AABB getBoundingBox(const MeshLib::Mesh& mesh);

    /**
     * Returns an array with the number of elements of each type in the given
     * mesh. On completion, n_element_types array contains the number of
     * elements of each of the seven supported types. The index to element type
     * conversion is this: 0: \#lines 1: \#triangles 2: \#quads 3: \#tetrahedra
     *        4: \#hexahedra
     *        5: \#pyramids
     *        6: \#prisms
     */
    static std::map<MeshElemType, unsigned> getNumberOfElementTypes(
        const MeshLib::Mesh& mesh);

    /// writes all numbers of element types
    static void writeAllNumbersOfElementTypes(const MeshLib::Mesh& mesh);

    /// writes out property vector information
    static void writePropertyVectorInformation(const MeshLib::Mesh& mesh);

    /// writes out mesh validation results
    /// Remark: MeshValidation can modify the original mesh
    static void writeMeshValidationResults(MeshLib::Mesh& mesh);
};

}  // namespace MeshLib
