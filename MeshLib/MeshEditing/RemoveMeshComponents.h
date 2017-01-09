/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <vector>

namespace MeshLib
{
class Mesh;
class Element;

/**
 * Removes mesh elements and returns a new mesh object. The original mesh is kept unchanged.
 * @param mesh                 an original mesh whose elements are removed
 * @param removed_element_ids  a vector of element indices to be removed
 * @param new_mesh_name        a new mesh name
 * @return a new mesh object
 */
MeshLib::Mesh* removeElements(const MeshLib::Mesh& mesh,
        const std::vector<std::size_t> &removed_element_ids, const std::string &new_mesh_name);

/**
 * Removes the mesh nodes (and connected elements) given in the nodes-list from the mesh.
 * @param mesh                 an original mesh whose elements are removed
 * @param removed_node_ids     a vector of node indices to be removed
 * @param new_mesh_name        a new mesh name
 * @return a new mesh object
 */
MeshLib::Mesh* removeNodes(const MeshLib::Mesh &mesh, const std::vector<std::size_t> &removed_node_ids, const std::string &new_mesh_name);

} // end namespace MeshLib
