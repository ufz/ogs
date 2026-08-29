// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#ifdef USE_PETSC
#include <petscsys.h>

#include "MeshLib/NodePartitionedMesh.h"
#endif

#include "MathLib/LinAlg/SparsityPattern.h"
#include "MeshLib/Elements/Utils.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"
#include "MeshToolsLib/MeshGenerators/QuadraticMeshGenerator.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/NumericsConfig.h"

// These tests build plain (non-partitioned) meshes with ComponentOrder::
// BY_COMPONENT and check the serial sparsity pattern. Under a PETSc build
// MeshComponentMap requires a NodePartitionedMesh and BY_LOCATION ordering, so
// the tests are DISABLED there. The expected nnzInRow values are the serial
// pattern's safe over-estimate (computeSparsityPatternNonPETSc counts each
// node's own DOFs both as "self" and again via the self-inclusive adjacency
// table).
#ifndef USE_PETSC
TEST(NumLib_SparsityPattern, SingleComponentLinearMesh)
#else
TEST(NumLib_SparsityPattern, DISABLED_SingleComponentLinearMesh)
#endif
{
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshToolsLib::MeshGenerator::generateLineMesh(3u, 1.));
    MeshLib::MeshSubset nodesSubset{*mesh, mesh->getNodes()};

    std::vector<MeshLib::MeshSubset> components{nodesSubset};
    NumLib::LocalToGlobalIndexMap dof_map(std::move(components),
                                          NumLib::ComponentOrder::BY_COMPONENT);

    GlobalSparsityPattern sp = NumLib::computeSparsityPattern(dof_map, *mesh);

    ASSERT_EQ(4, sp.numberOfRows());
    EXPECT_EQ(3, sp.nnzInRow(0));
    EXPECT_EQ(4, sp.nnzInRow(1));
    EXPECT_EQ(4, sp.nnzInRow(2));
    EXPECT_EQ(3, sp.nnzInRow(3));
}

#ifndef USE_PETSC
TEST(NumLib_SparsityPattern, SingleComponentQuadraticMesh)
#else
TEST(NumLib_SparsityPattern, DISABLED_SingleComponentQuadraticMesh)
#endif
{
    std::unique_ptr<MeshLib::Mesh> linear_mesh(
        MeshToolsLib::MeshGenerator::generateLineMesh(3u, 1.));
    std::unique_ptr<MeshLib::Mesh> mesh(MeshToolsLib::createQuadraticOrderMesh(
        *linear_mesh, false /* add centre node */));
    MeshLib::MeshSubset nodesSubset{*mesh, mesh->getNodes()};

    std::vector<MeshLib::MeshSubset> components{nodesSubset};
    NumLib::LocalToGlobalIndexMap dof_map(std::move(components),
                                          NumLib::ComponentOrder::BY_COMPONENT);

    GlobalSparsityPattern sp = NumLib::computeSparsityPattern(dof_map, *mesh);

    ASSERT_EQ(7, sp.numberOfRows());
    EXPECT_EQ(4, sp.nnzInRow(0));
    EXPECT_EQ(6, sp.nnzInRow(1));
    EXPECT_EQ(6, sp.nnzInRow(2));
    EXPECT_EQ(4, sp.nnzInRow(3));
    EXPECT_EQ(4, sp.nnzInRow(4));
    EXPECT_EQ(4, sp.nnzInRow(5));
    EXPECT_EQ(4, sp.nnzInRow(6));
}

#ifndef USE_PETSC
TEST(NumLib_SparsityPattern, MultipleComponentsLinearMesh)
#else
TEST(NumLib_SparsityPattern, DISABLED_MultipleComponentsLinearMesh)
#endif
{
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshToolsLib::MeshGenerator::generateLineMesh(3u, 1.));
    MeshLib::MeshSubset nodesSubset{*mesh, mesh->getNodes()};

    std::vector<MeshLib::MeshSubset> components{nodesSubset, nodesSubset};
    NumLib::LocalToGlobalIndexMap dof_map(std::move(components),
                                          NumLib::ComponentOrder::BY_COMPONENT);

    GlobalSparsityPattern sp = NumLib::computeSparsityPattern(dof_map, *mesh);

    ASSERT_EQ(8, sp.numberOfRows());
    auto const n = static_cast<GlobalIndexType>(mesh->getNumberOfNodes());
    for (GlobalIndexType i = 0; i < 2; i++)
    {
        EXPECT_EQ(6, sp.nnzInRow(i * n + 0));
        EXPECT_EQ(8, sp.nnzInRow(i * n + 1));
        EXPECT_EQ(8, sp.nnzInRow(i * n + 2));
        EXPECT_EQ(6, sp.nnzInRow(i * n + 3));
    }
}

#ifndef USE_PETSC
TEST(NumLib_SparsityPattern, MultipleComponentsLinearQuadraticMesh)
#else
TEST(NumLib_SparsityPattern, DISABLED_MultipleComponentsLinearQuadraticMesh)
#endif
{
    std::unique_ptr<MeshLib::Mesh> linear_mesh(
        MeshToolsLib::MeshGenerator::generateLineMesh(3u, 1.));
    std::unique_ptr<MeshLib::Mesh> mesh(MeshToolsLib::createQuadraticOrderMesh(
        *linear_mesh, false /* add centre node */));
    auto base_nodes = MeshLib::getBaseNodes(mesh->getElements());
    auto baseNodesSubset =
        std::make_unique<MeshLib::MeshSubset const>(*mesh, base_nodes);
    auto allNodesSubset =
        std::make_unique<MeshLib::MeshSubset const>(*mesh, mesh->getNodes());

    std::vector<MeshLib::MeshSubset> components{*baseNodesSubset,
                                                *allNodesSubset};
    NumLib::LocalToGlobalIndexMap dof_map(std::move(components),
                                          NumLib::ComponentOrder::BY_COMPONENT);

    GlobalSparsityPattern sp = NumLib::computeSparsityPattern(dof_map, *mesh);

    ASSERT_EQ(11, sp.numberOfRows());
    // 1st component
    EXPECT_EQ(7, sp.nnzInRow(0));
    EXPECT_EQ(10, sp.nnzInRow(1));
    EXPECT_EQ(10, sp.nnzInRow(2));
    EXPECT_EQ(7, sp.nnzInRow(3));
    // 2nd component
    EXPECT_EQ(7, sp.nnzInRow(4));
    EXPECT_EQ(10, sp.nnzInRow(5));
    EXPECT_EQ(10, sp.nnzInRow(6));
    EXPECT_EQ(7, sp.nnzInRow(7));
    EXPECT_EQ(6, sp.nnzInRow(8));
    EXPECT_EQ(6, sp.nnzInRow(9));
    EXPECT_EQ(6, sp.nnzInRow(10));
}

#ifdef USE_PETSC
// PETSc-native check of computeSparsityPatternPETSc. In contrast to the serial
// path above, it produces the exact, deduplicated CSR pattern from mesh
// adjacency (each nonzero counted once). Registered without an MPI_ name, so it
// runs under the single-rank testrunner invocation; a single-partition
// NodePartitionedMesh is only valid there. The multi-rank / ghost-node path is
// covered by the MPI matrix-interface tests.
TEST(NumLib_SparsityPattern, PETScSingleComponentLinearMesh)
{
    int comm_size;
    MPI_Comm_size(PETSC_COMM_WORLD, &comm_size);
    if (comm_size != 1)
    {
        GTEST_SKIP() << "A single-partition mesh is only valid on one rank.";
    }

    std::unique_ptr<MeshLib::Mesh> linear_mesh(
        MeshToolsLib::MeshGenerator::generateLineMesh(3u, 1.));
    MeshLib::NodePartitionedMesh const mesh{*linear_mesh};

    MeshLib::MeshSubset nodesSubset{mesh, mesh.getNodes()};
    std::vector<MeshLib::MeshSubset> components{nodesSubset};
    NumLib::LocalToGlobalIndexMap dof_map(std::move(components),
                                          NumLib::ComponentOrder::BY_LOCATION);

    GlobalSparsityPattern const sp =
        NumLib::computeSparsityPattern(dof_map, mesh);

    // Exact adjacency pattern (deduplicated), not the serial over-estimate:
    // a 3-cell line has nonzero counts 2, 3, 3, 2.
    ASSERT_EQ(4, sp.numberOfRows());
    EXPECT_EQ(2, sp.nnzInRow(0));
    EXPECT_EQ(3, sp.nnzInRow(1));
    EXPECT_EQ(3, sp.nnzInRow(2));
    EXPECT_EQ(2, sp.nnzInRow(3));

    EXPECT_EQ((std::vector<GlobalIndexType>{0, 2, 5, 8, 10}), sp.row_ptr);
    EXPECT_EQ((std::vector<GlobalIndexType>{0, 1, 0, 1, 2, 1, 2, 3, 2, 3}),
              sp.col_idx);
}
#endif
