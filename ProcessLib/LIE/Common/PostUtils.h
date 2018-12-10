/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <map>
#include <memory>
#include <vector>

#include "MeshLib/Mesh.h"
namespace MeshLib
{
class Node;
}

namespace ProcessLib
{
namespace LIE
{
/// A tool for post-processing results from the LIE approach
///
/// The tool creates a new mesh containing duplicated fracture nodes
/// to represent geometric discontinuities in visualization.
class PostProcessTool final
{
public:
    PostProcessTool(
        MeshLib::Mesh const& org_mesh,
        std::vector<int> const& vec_fracture_mat_IDs,
        std::vector<std::vector<MeshLib::Node*>> const& vec_vec_fracture_nodes,
        std::vector<std::vector<MeshLib::Element*>> const&
            vec_vec_fracture_matrix_elements,
        std::vector<std::pair<std::size_t, std::vector<int>>> const&
            vec_branch_nodeID_matIDs,
        std::vector<std::pair<std::size_t, std::vector<int>>> const&
            vec_junction_nodeID_matIDs);

    MeshLib::Mesh const& getOutputMesh() const { return *_output_mesh; }

private:
    template <typename T>
    void createProperties();
    template <typename T>
    void copyProperties();
    void calculateTotalDisplacement(unsigned const n_fractures,
                                    unsigned const n_junctions);

    MeshLib::Mesh const& _org_mesh;
    std::unique_ptr<MeshLib::Mesh> _output_mesh;
    std::map<std::size_t, std::vector<std::size_t>> _map_dup_newNodeIDs;
};

}  // namespace LIE
}  // namespace ProcessLib
