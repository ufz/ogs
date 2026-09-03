// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ComputeSparsityPattern.h"

#include <numeric>
#include <range/v3/algorithm/sort.hpp>
#include <range/v3/algorithm/unique.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>

#include "LocalToGlobalIndexMap.h"
#include "MeshLib/NodeAdjacencyTable.h"

#ifdef USE_PETSC
#include "MeshLib/NodePartitionedMesh.h"

GlobalSparsityPattern computeSparsityPatternPETSc(
    NumLib::LocalToGlobalIndexMap const& dof_table, MeshLib::Mesh const& mesh)
{
    assert(dynamic_cast<MeshLib::NodePartitionedMesh const*>(&mesh));

    MeshLib::NodeAdjacencyTable const node_adjacency_table(mesh);

    // getGlobalIndices() returns raw PETSc global indices:
    //   idx >= 0  — owned by this rank, actual PETSc global row/col index
    //   idx < 0   — ghost node on this rank (owned by another rank);
    //               the actual global index is stored negated
    auto const global_idcs =
        MeshLib::views::meshLocations(mesh, MeshLib::MeshItemType::Node) |
        ranges::views::transform([&](auto&& l)
                                 { return dof_table.getGlobalIndices(l); }) |
        ranges::to<std::vector>();

    auto const n_local =
        static_cast<GlobalIndexType>(dof_table.dofSizeWithoutGhosts());

    // Each rank owns a contiguous range of global indices; global_start is the
    // first index owned here, obtained by an exclusive prefix-sum of the local
    // sizes. It maps a global row index to the corresponding local row index.
    //
    // MPIU_INT is not an MPI standard datatype but a PETSc macro (defined in
    // petscsys.h) that expands to the MPI datatype matching the width of
    // PetscInt, i.e. MPI_INT or MPI_INT64_T depending on the PETSc
    // configuration. It is used rather than the BaseLib::MPI wrappers because
    // BaseLib::MPI::mpiType() has no mapping for a 64-bit PetscInt and would
    // silently pick the wrong datatype.
    GlobalIndexType global_start = 0;
    MPI_Exscan(&n_local, &global_start, 1, MPIU_INT, MPI_SUM, PETSC_COMM_WORLD);
    // Total global DOF count needed to decode the special ghost encoding where
    // global index 0 is stored as -num_global_dof.
    GlobalIndexType num_global_dof = 0;
    MPI_Allreduce(&n_local, &num_global_dof, 1, MPIU_INT, MPI_SUM,
                  PETSC_COMM_WORLD);

    // Collect all (row, col) pairs into per-local-row lists of global column
    // indices; they are sorted and deduplicated below.
    //
    // A row holds only a few dozen columns, so a flat vector plus one
    // sort/unique per row is used instead of std::set: appending is a
    // contiguous write without a per-element node allocation, and the sort
    // is on data that already sits in cache. The result is the same -- sorted
    // and unique -- as std::set gives, which is what the CSR pattern needs.
    std::vector<std::vector<GlobalIndexType>> col_lists(n_local);

    // Decode a raw PETSc global index (possibly negative ghost encoding) to the
    // actual non-negative global index.
    auto const decode_col = [&](GlobalIndexType const raw) -> GlobalIndexType
    {
        if (raw >= 0)
        {
            return raw;
        }
        // Ghost encoding: actual index k is stored as -k, except k==0 which is
        // stored as -num_global_dof.
        return (-raw == num_global_dof) ? 0 : -raw;
    };

    // Collect a (row, col) pair into the per-row column set.
    //   row >= 0 (non-negative)  — owned by this rank
    //   row < 0                  — ghost row, not owned here → skip
    auto const collect =
        [&](GlobalIndexType const row, GlobalIndexType const col)
    {
        if (row < 0)
        {
            return;
        }
        // A non-negative raw index is owned by this rank, and an owned index
        // lies in this rank's contiguous range [global_start, global_start +
        // n_local), so the local row is always in range here.
        GlobalIndexType const local_row = row - global_start;
        assert(local_row >= 0 && local_row < n_local);
        col_lists[local_row].push_back(decode_col(col));
    };

    // Standard mesh-adjacency contribution.
    // Note: getAdjacentNodes() already includes node n itself (the adjacency
    // table is built from all element nodes, including self), so there is no
    // need for a separate self-node loop.
    for (std::size_t n = 0; n < mesh.getNumberOfNodes(); ++n)
    {
        auto const& adj_nodes = node_adjacency_table.getAdjacentNodes(n);
        for (auto const row_idx : global_idcs[n])
        {
            for (auto const adj : adj_nodes)
            {
                for (auto const col_idx : global_idcs[adj])
                {
                    collect(row_idx, col_idx);
                }
            }
        }
    }

    // Sort and deduplicate each row, so that col_idx ends up sorted within each
    // row as PETScSparsityPattern documents.
    std::size_t number_of_nonzeros = 0;
    for (auto& cols : col_lists)
    {
        ranges::sort(cols);
        cols.erase(ranges::unique(cols), cols.end());
        number_of_nonzeros += cols.size();
    }

    // Build the CSR sparsity pattern from the collected column lists.
    GlobalSparsityPattern sparsity_pattern;
    sparsity_pattern.row_ptr.resize(n_local + 1);
    sparsity_pattern.row_ptr[0] = 0;
    sparsity_pattern.col_idx.reserve(number_of_nonzeros);

    for (GlobalIndexType local_row = 0; local_row < n_local; ++local_row)
    {
        auto const& cols = col_lists[local_row];
        sparsity_pattern.col_idx.insert(sparsity_pattern.col_idx.end(),
                                        cols.begin(), cols.end());
        sparsity_pattern.row_ptr[local_row + 1] =
            sparsity_pattern.row_ptr[local_row] +
            static_cast<GlobalIndexType>(cols.size());
    }

    return sparsity_pattern;
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

    GlobalSparsityPattern sparsity_pattern;
    sparsity_pattern.number_non_zeros_per_row.assign(
        dof_table.dofSizeWithGhosts(), 0);

    // Map adjacent mesh nodes to "adjacent global indices".
    for (std::size_t n = 0; n < mesh.getNumberOfNodes(); ++n)
    {
        auto const& an = node_adjacency_table.getAdjacentNodes(n);
        auto const n_self_dof = global_idcs[n].size();
        auto const n_connected_dof = std::accumulate(
            cbegin(an), cend(an), 0, [&](auto const result, auto const i)
            { return result + global_idcs[i].size(); });
        auto const n_dof = n_self_dof + n_connected_dof;
        for (auto global_index : global_idcs[n])
        {
            sparsity_pattern.number_non_zeros_per_row[global_index] = n_dof;
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
