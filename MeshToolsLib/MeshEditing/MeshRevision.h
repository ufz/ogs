/**
 * \file
 * \author Karsten Rink
 * \date   2014-02-14
 * \brief  Definition of the MeshRevision class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>
#include <limits>
#include <string>
#include <vector>

// forward declaration
namespace MeshLib
{
class Mesh;
class Node;
}  // namespace MeshLib

namespace MeshToolsLib
{
/**
 * Collapses nodes with a distance smaller min_distance and
 * reduces elements accordingly.
 */
class MeshRevision
{
public:
    /**
     * Constructor
     * @param mesh The mesh which is being revised. Note that node IDs in mesh
     * are changed during computation but are reset after the algorithms
     * implemented here are finished
     */
    explicit MeshRevision(MeshLib::Mesh& mesh);

    /// Returns the number of potentially collapsible nodes
    unsigned getNumberOfCollapsibleNodes(
        double eps = std::numeric_limits<double>::epsilon()) const;

    /// Designates nodes to be collapsed by setting their ID to the index of the
    /// node they will get merged with.
    std::vector<std::size_t> collapseNodeIndices(double eps) const;

    /**
     * Create a new mesh where all nodes with a distance < eps from each other
     * are collapsed.
     * Elements are adjusted accordingly and elements with nonplanar faces are
     * subdivided into geometrically correct elements.
     * @param new_mesh_name New name.
     * @param eps           Minimum distance for nodes not to be collapsed
     * @param min_elem_dim  Minimum dimension of elements to be inserted into
     *                      new mesh (i.e. min_elem_dim=3 will prevent the new
     *                      mesh to contain 2D elements)
     */
    MeshLib::Mesh* simplifyMesh(const std::string& new_mesh_name, double eps,
                                unsigned min_elem_dim = 1) const;

private:
    /// Constructs a new node vector for the resulting mesh by removing all
    /// nodes whose ID indicates they need to be merged/removed.
    std::vector<MeshLib::Node*> constructNewNodesArray(
        const std::vector<std::size_t>& id_map) const;

    /// The original mesh used for constructing the class
    MeshLib::Mesh& _mesh;
};

}  // namespace MeshToolsLib
