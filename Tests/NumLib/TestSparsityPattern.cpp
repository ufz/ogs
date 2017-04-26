/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshGenerators/QuadraticMeshGenerator.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/NumericsConfig.h"


#ifndef USE_PETSC
TEST(NumLib_SparsityPattern, SingleComponentLinearMesh)
#else
TEST(NumLib_SparsityPattern, DISABLED_SingleComponentLinearMesh)
#endif
{
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::MeshGenerator::generateLineMesh(3u, 1.));
    std::unique_ptr<MeshLib::MeshSubset const> nodesSubset(
        new MeshLib::MeshSubset(*mesh, &mesh->getNodes()));

    std::vector<MeshLib::MeshSubsets> components;
    components.emplace_back(nodesSubset.get());
    NumLib::LocalToGlobalIndexMap dof_map(
                      std::move(components),
                      NumLib::ComponentOrder::BY_COMPONENT);

    GlobalSparsityPattern sp = NumLib::computeSparsityPattern(dof_map, *mesh.get());

    ASSERT_EQ(4u, sp.size());
    EXPECT_EQ(2u, sp[0]);
    EXPECT_EQ(3u, sp[1]);
    EXPECT_EQ(3u, sp[2]);
    EXPECT_EQ(2u, sp[3]);
}


#ifndef USE_PETSC
TEST(NumLib_SparsityPattern, SingleComponentQuadraticMesh)
#else
TEST(NumLib_SparsityPattern, DISABLED_SingleComponentQuadraticMesh)
#endif
{
    std::unique_ptr<MeshLib::Mesh> linear_mesh(
        MeshLib::MeshGenerator::generateLineMesh(3u, 1.));
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::createQuadraticOrderMesh(*linear_mesh));
    std::unique_ptr<MeshLib::MeshSubset const> nodesSubset(
        new MeshLib::MeshSubset(*mesh, &mesh->getNodes()));

    std::vector<MeshLib::MeshSubsets> components;
    components.emplace_back(MeshLib::MeshSubsets{nodesSubset.get()});
    NumLib::LocalToGlobalIndexMap dof_map(
                      std::move(components),
                      NumLib::ComponentOrder::BY_COMPONENT);

    GlobalSparsityPattern sp = NumLib::computeSparsityPattern(dof_map, *mesh.get());

    ASSERT_EQ(7u, sp.size());
    EXPECT_EQ(3u, sp[0]);
    EXPECT_EQ(5u, sp[1]);
    EXPECT_EQ(5u, sp[2]);
    EXPECT_EQ(3u, sp[3]);
    EXPECT_EQ(3u, sp[4]);
    EXPECT_EQ(3u, sp[5]);
    EXPECT_EQ(3u, sp[6]);
}


#ifndef USE_PETSC
TEST(NumLib_SparsityPattern, MultipleComponentsLinearMesh)
#else
TEST(NumLib_SparsityPattern, DISABLED_MultipleComponentsLinearMesh)
#endif
{
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::MeshGenerator::generateLineMesh(3u, 1.));
    std::unique_ptr<MeshLib::MeshSubset const> nodesSubset(
        new MeshLib::MeshSubset(*mesh, &mesh->getNodes()));

    std::vector<MeshLib::MeshSubsets> components;
    components.emplace_back(nodesSubset.get());
    components.emplace_back(nodesSubset.get());
    NumLib::LocalToGlobalIndexMap dof_map(
                      std::move(components),
                      NumLib::ComponentOrder::BY_COMPONENT);

    GlobalSparsityPattern sp = NumLib::computeSparsityPattern(dof_map, *mesh.get());

    ASSERT_EQ(8u, sp.size());
    for (int i=0; i<2; i++)
    {
        EXPECT_EQ(4u, sp[i*mesh->getNumberOfNodes() + 0]);
        EXPECT_EQ(6u, sp[i*mesh->getNumberOfNodes() + 1]);
        EXPECT_EQ(6u, sp[i*mesh->getNumberOfNodes() + 2]);
        EXPECT_EQ(4u, sp[i*mesh->getNumberOfNodes() + 3]);
    }
}


#ifndef USE_PETSC
TEST(NumLib_SparsityPattern, MultipleComponentsLinearQuadraticMesh)
#else
TEST(NumLib_SparsityPattern, DISABLED_MultipleComponentsLinearQuadraticMesh)
#endif
{
    std::unique_ptr<MeshLib::Mesh> linear_mesh(
        MeshLib::MeshGenerator::generateLineMesh(3u, 1.));
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::createQuadraticOrderMesh(*linear_mesh));
    auto base_nodes = MeshLib::getBaseNodes(mesh->getElements());
    std::unique_ptr<MeshLib::MeshSubset const> baseNodesSubset(
        new MeshLib::MeshSubset(*mesh, &base_nodes));
    std::unique_ptr<MeshLib::MeshSubset const> allNodesSubset(
        new MeshLib::MeshSubset(*mesh, &mesh->getNodes()));

    std::vector<MeshLib::MeshSubsets> components;
    components.emplace_back(baseNodesSubset.get());
    components.emplace_back(allNodesSubset.get());
    NumLib::LocalToGlobalIndexMap dof_map(
                      std::move(components),
                      NumLib::ComponentOrder::BY_COMPONENT);

    GlobalSparsityPattern sp = NumLib::computeSparsityPattern(dof_map, *mesh.get());

    ASSERT_EQ(11u, sp.size());
    // 1st component
    EXPECT_EQ(5u, sp[0]);
    EXPECT_EQ(8u, sp[1]);
    EXPECT_EQ(8u, sp[2]);
    EXPECT_EQ(5u, sp[3]);
    // 2nd component
    EXPECT_EQ(5u, sp[4]);
    EXPECT_EQ(8u, sp[5]);
    EXPECT_EQ(8u, sp[6]);
    EXPECT_EQ(5u, sp[7]);
    EXPECT_EQ(5u, sp[8]);
    EXPECT_EQ(5u, sp[9]);
    EXPECT_EQ(5u, sp[10]);
}

