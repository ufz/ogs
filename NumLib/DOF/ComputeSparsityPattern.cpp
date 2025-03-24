/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ComputeSparsityPattern.h"

#include <numeric>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>

#include "LocalToGlobalIndexMap.h"
#include "MeshLib/Location.h"
#include "MeshLib/NodeAdjacencyTable.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"

GlobalSparsityPattern computeSparsityPatternPETSc(
    NumLib::LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh)
{
    assert(dynamic_cast<MeshLib::NodePartitionedMesh const*>(&mesh));
    auto const& npmesh =
        *static_cast<MeshLib::NodePartitionedMesh const*>(&mesh);

    auto const max_nonzeroes = dof_table.getNumberOfGlobalComponents() *
                               npmesh.getMaximumNConnectedNodesToNode();

    // The sparsity pattern is misused here in the sense that it will only
    // contain a single value.
    return GlobalSparsityPattern(1, max_nonzeroes);
}
#else
GlobalSparsityPattern computeSparsityPatternNonPETSc(
    NumLib::LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh)
{
    MeshLib::NodeAdjacencyTable const node_adjacency_table(mesh);

    // A mapping   mesh node id -> global indices
    // It acts as a cache for dof table queries.
    auto const global_idcs =
        MeshLib::views::meshLocations(mesh, MeshLib::MeshItemType::Node) |
        ranges::views::transform([&](auto&& l)
                                 { return dof_table.getGlobalIndices(l); }) |
        ranges::to<std::vector>();

    GlobalSparsityPattern sparsity_pattern(dof_table.dofSizeWithGhosts());

    // Map adjacent mesh nodes to "adjacent global indices".
    for (std::size_t n = 0; n < mesh.getNumberOfNodes(); ++n)
    {
        auto const& an = node_adjacency_table.getAdjacentNodes(n);
        auto const n_self_dof = global_idcs[n].size();
        auto const n_connected_dof =
            std::accumulate(cbegin(an), cend(an), 0,
                            [&](auto const result, auto const i)
                            { return result + global_idcs[i].size(); });
        auto const n_dof = n_self_dof + n_connected_dof;
        for (auto global_index : global_idcs[n])
        {
            sparsity_pattern[global_index] = n_dof;
        }
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

}  // namespace NumLib
