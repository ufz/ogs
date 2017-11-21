/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <map>
#include <memory>

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace LIE
{
/// A tool for post-processing results from the LIE approach
///
/// The tool creates a new mesh containing duplicated fracture nodes
/// to represent geometric discontinuities in visualization.
class PostProcessTool
{
public:
    PostProcessTool(
        MeshLib::Mesh const& org_mesh,
        std::vector<std::vector<MeshLib::Node*>> const& vec_vec_fracture_nodes,
        std::vector<std::vector<MeshLib::Element*>> const&
            vec_vec_fracutre_matrix_elements);

    MeshLib::Mesh const& getOutputMesh() const { return *_output_mesh; }

private:
    template <typename T>
    void createProperties();
    template <typename T>
    void copyProperties();
    void calculateTotalDisplacement(unsigned const n_fractures);

    MeshLib::Mesh const& _org_mesh;
    std::unique_ptr<MeshLib::Mesh> _output_mesh;
    std::map<std::size_t, std::size_t> _map_dup_newNodeIDs;
};

}  // namespace LIE
}  // namespace ProcessLib
