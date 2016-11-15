/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ComputeSparsityPattern.h"

#include "LocalToGlobalIndexMap.h"
#include "MeshLib/NodeAdjacencyTable.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"

GlobalSparsityPattern computeSparsityPatternPETSc(
    NumLib::LocalToGlobalIndexMap const& dof_table,
    MeshLib::Mesh const& mesh)
{
    assert(dynamic_cast<MeshLib::NodePartitionedMesh const*>(&mesh));
    auto const& npmesh =
        *static_cast<MeshLib::NodePartitionedMesh const*>(&mesh);

    auto const max_nonzeroes =   dof_table.getNumberOfComponents()
                               * npmesh.getMaximumNConnectedNodesToNode();

    // The sparsity pattern is misused here in the sense that it will only
    // contain a single value.
    return GlobalSparsityPattern(1, max_nonzeroes);
}
#else
GlobalSparsityPattern computeSparsityPatternNonPETSc(
    NumLib::LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh)
{
    MeshLib::NodeAdjacencyTable node_adjacency_table;
    node_adjacency_table.createTable(mesh.getNodes());

    // A mapping   mesh node id -> global indices
    // It acts as a cache for dof table queries.
    std::vector<std::vector<GlobalIndexType>> global_idcs;

    global_idcs.reserve(mesh.getNumberOfNodes());
    for (std::size_t n = 0; n < mesh.getNumberOfNodes(); ++n)
    {
        MeshLib::Location l(mesh.getID(), MeshLib::MeshItemType::Node, n);
        global_idcs.push_back(dof_table.getGlobalIndices(l));
    }

    GlobalSparsityPattern sparsity_pattern(dof_table.dofSizeWithGhosts());

    // Map adjacent mesh nodes to "adjacent global indices".
    for (std::size_t n = 0; n < mesh.getNumberOfNodes(); ++n)
    {
        unsigned n_connected_dof = 0;
        for (auto an : node_adjacency_table.getAdjacentNodes(n))
            n_connected_dof += global_idcs[an].size();
        for (auto global_index : global_idcs[n])
            sparsity_pattern[global_index] = n_connected_dof;
    }

    return sparsity_pattern;
}
#endif

namespace NumLib
{
GlobalSparsityPattern computeSparsityPattern(
    LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh)
{
#ifdef USE_PETSC
    return computeSparsityPatternPETSc(dof_table, mesh);
#else
    return computeSparsityPatternNonPETSc(dof_table, mesh);
#endif
}

}
