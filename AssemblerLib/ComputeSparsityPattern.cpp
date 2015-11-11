/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ComputeSparsityPattern.h"

namespace AssemblerLib
{

SparsityPattern
computeSparsityPattern(
        MeshLib::NodeAdjacencyTable const& node_adjacency_table,
        LocalToGlobalIndexMap const& dof_table,
        MeshLib::Mesh const& mesh
        )
{
    std::vector<std::vector<GlobalIndexType> > global_idcs;
    global_idcs.reserve(mesh.getNNodes());
    for (std::size_t n=0; n<mesh.getNNodes(); ++n)
    {
        MeshLib::Location l(mesh.getID(), MeshLib::MeshItemType::Node, n);
        global_idcs.push_back(dof_table.getGlobalIndices(l));
    }

    SparsityPattern sparsity_pattern(dof_table.dofSize());

    for (std::size_t n=0; n<mesh.getNNodes(); ++n)
    {
        auto const& node_ids = node_adjacency_table.getAdjacentNodes(n);
        for (auto an : node_ids) {
            auto const& row_ids = global_idcs[an];
            for (auto r : row_ids) {
                sparsity_pattern[r] += row_ids.size();
            }
        }
    }

    return sparsity_pattern;
}

}
