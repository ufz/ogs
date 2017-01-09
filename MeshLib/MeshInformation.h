/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-28
 * \brief  Definition of the MeshInformation class
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>
#include <string>
#include <limits>

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
    /// Returns the smallest and largest value of a scalar array with the given name.
    template<typename T>
    static std::pair<T, T> const
        getValueBounds(MeshLib::Mesh const& mesh, std::string const& name)
    {
        auto const* const data_vec =
            mesh.getProperties().getPropertyVector<T>(name);
        if (!data_vec)
            return {std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};
        if (data_vec->empty()) {
            INFO("Mesh does not contain values for the property \"%s\".", name.c_str());
            return {std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};
        }
        auto vec_bounds = std::minmax_element(data_vec->cbegin(), data_vec->cend());
        return {*(vec_bounds.first), *(vec_bounds.second)};
    }

    /// Returns the bounding box of the mesh.
    static const GeoLib::AABB getBoundingBox(const MeshLib::Mesh &mesh);

    /**
     * Returns an array with the number of elements of each type in the given mesh.
     * On completion, n_element_types array contains the number of elements of each of the seven
     * supported types. The index to element type conversion is this:
     *        0: \#lines
     *        1: \#triangles
     *        2: \#quads
     *        3: \#tetrahedra
     *        4: \#hexahedra
     *        5: \#pyramids
     *        6: \#prisms
     */
    static const std::array<unsigned, 7> getNumberOfElementTypes(const MeshLib::Mesh &mesh);


};

}
