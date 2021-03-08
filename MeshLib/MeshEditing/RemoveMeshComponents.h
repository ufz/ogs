/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
class Node;

/**
 * Removes mesh elements and returns a new mesh object. The original mesh is kept unchanged.
 * @param mesh                 an original mesh whose elements are removed
 * @param removed_element_ids  a vector of element indices to be removed
 * @param new_mesh_name        a new mesh name
 * @return a new mesh object
 */
MeshLib::Mesh* removeElements(
    const MeshLib::Mesh& mesh,
    const std::vector<std::size_t>& removed_element_ids,
    const std::string& new_mesh_name);

/**
 * Removes the mesh nodes (and connected elements) given in the nodes-list from
 * the mesh.
 * @param mesh                 an original mesh whose elements are removed
 * @param del_nodes_idx        a vector of node indices to be removed
 * @param new_mesh_name        a new mesh name
 * @return a new mesh object
 */
MeshLib::Mesh* removeNodes(const MeshLib::Mesh& mesh,
                           const std::vector<std::size_t>& del_nodes_idx,
                           const std::string& new_mesh_name);

/// Marks nodes not used by any of the elements.
std::vector<bool> markUnusedNodes(std::vector<Element*> const& elements,
                                  std::vector<Node*> const& nodes);

/// Deallocates and removes nodes marked true.
void removeMarkedNodes(std::vector<bool> const& nodes_to_delete,
                       std::vector<Node*>& nodes);
} // end namespace MeshLib
