// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#ifdef USE_PETSC
#include "BaseLib/MPI.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/MeshSubset.h"
#include "MeshLib/Node.h"
#include "MeshLib/NodePartitionedMesh.h"
#include "NumLib/DOF/MeshComponentMap.h"

namespace
{
/// A line mesh of 6 nodes at positions 0..5 on 3 ranks: rank r owns the nodes
/// 2r and 2r+1 and has the neighbouring nodes as ghosts. The global node ids
/// are the positions, or reversed, which breaks the partition block layout as
/// in subdomain meshes.
std::unique_ptr<MeshLib::NodePartitionedMesh> createPartitionedLineMesh(
    int const rank, bool const reversed_global_ids)
{
    std::vector<std::size_t> positions{2u * rank, 2u * rank + 1};
    if (rank > 0)
    {
        positions.push_back(2 * rank - 1);
    }
    if (rank < 2)
    {
        positions.push_back(2 * rank + 2);
    }

    std::vector<MeshLib::Node*> nodes;
    std::vector<std::size_t> global_ids;
    for (std::size_t i = 0; i < positions.size(); ++i)
    {
        nodes.push_back(new MeshLib::Node(positions[i], 0, 0, i));
        global_ids.push_back(reversed_global_ids ? 5 - positions[i]
                                                 : positions[i]);
    }
    std::vector<MeshLib::Element*> elements;
    for (std::size_t a = 0; a < nodes.size(); ++a)
    {
        for (std::size_t b = 0; b < nodes.size(); ++b)
        {
            if (positions[b] == positions[a] + 1)
            {
                elements.push_back(new MeshLib::Line(
                    std::array<MeshLib::Node*, 2>{nodes[a], nodes[b]}));
            }
        }
    }

    return std::make_unique<MeshLib::NodePartitionedMesh>(
        "line", nodes, global_ids, elements, MeshLib::Properties{}, 6, 6, 2,
        std::vector<std::size_t>{0, 2, 4, 6},
        std::vector<std::size_t>{0, 0, 0, 0});
}

/// The d.o.f. of a node is numbered by its owning rank in node order, so its
/// global index is the node's position, on the owning rank and as a ghost.
void checkGlobalIndicesArePositions(bool const reversed_global_ids)
{
    BaseLib::MPI::Mpi const mpi;
    if (mpi.size != 3)
    {
        GTEST_SKIP() << "Test requires exactly 3 MPI ranks";
    }

    auto const mesh = createPartitionedLineMesh(mpi.rank, reversed_global_ids);
    std::vector<MeshLib::MeshSubset> const components{
        {*mesh, mesh->getNodes()}};
    NumLib::MeshComponentMap const map(components,
                                       NumLib::ComponentOrder::BY_LOCATION);

    ASSERT_EQ(2, map.dofSizeWithoutGhosts());
    ASSERT_EQ(mesh->getNumberOfNodes(), map.dofSizeWithGhosts());
    for (std::size_t i = 0; i < 2; ++i)
    {
        EXPECT_EQ((*mesh->getNode(i))[0],
                  map.getGlobalIndex(
                      {mesh->getID(), MeshLib::MeshItemType::Node, i}, 0));
    }
    std::vector<GlobalIndexType> expected_ghost_indices;
    for (std::size_t i = 2; i < mesh->getNumberOfNodes(); ++i)
    {
        expected_ghost_indices.push_back((*mesh->getNode(i))[0]);
    }
    EXPECT_THAT(map.getGhostIndices(),
                ::testing::ElementsAreArray(expected_ghost_indices));
}
}  // namespace

TEST(MPI_MeshComponentMapParallel, GlobalIdsInPartitionBlockLayout)
{
    checkGlobalIndicesArePositions(false);
}

TEST(MPI_MeshComponentMapParallel, GlobalIdsNotInPartitionBlockLayout)
{
    checkGlobalIndicesArePositions(true);
}

#endif  // USE_PETSC
