/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <gtest/gtest.h>

#include "MeshLib/Elements/Utils.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshGenerators/QuadraticeMeshGenerator.h"

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/DOF/ComputeSparsityPattern.h"
#include "NumLib/NumericsConfig.h"


TEST(NumLib_SparsityPattern, SingleComponentLinearMesh)
{
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::MeshGenerator::generateLineMesh(3u, 1.));
    std::unique_ptr<MeshLib::MeshSubset const> nodesSubset(
        new MeshLib::MeshSubset(*mesh, &mesh->getNodes()));

    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> components;
    components.emplace_back(new MeshLib::MeshSubsets{nodesSubset.get()});
    NumLib::LocalToGlobalIndexMap dof_map(
                      std::move(components),
                      NumLib::ComponentOrder::BY_COMPONENT);

    GlobalSparsityPattern sp = NumLib::computeSparsityPattern(dof_map, *mesh.get());

    ASSERT_EQ(4u, sp.size());
    ASSERT_EQ(2u, sp[0]);
    ASSERT_EQ(3u, sp[1]);
    ASSERT_EQ(3u, sp[2]);
    ASSERT_EQ(2u, sp[3]);
}


TEST(NumLib_SparsityPattern, SingleComponentQuadraticMesh)
{
    std::unique_ptr<MeshLib::Mesh> linear_mesh(
        MeshLib::MeshGenerator::generateLineMesh(3u, 1.));
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::createQuadraticOrderMesh(*linear_mesh));
    std::unique_ptr<MeshLib::MeshSubset const> nodesSubset(
        new MeshLib::MeshSubset(*mesh, &mesh->getNodes()));

    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> components;
    components.emplace_back(new MeshLib::MeshSubsets{nodesSubset.get()});
    NumLib::LocalToGlobalIndexMap dof_map(
                      std::move(components),
                      NumLib::ComponentOrder::BY_COMPONENT);

    GlobalSparsityPattern sp = NumLib::computeSparsityPattern(dof_map, *mesh.get());

    ASSERT_EQ(7u, sp.size());
    ASSERT_EQ(3u, sp[0]);
    ASSERT_EQ(5u, sp[1]);
    ASSERT_EQ(5u, sp[2]);
    ASSERT_EQ(3u, sp[3]);
    ASSERT_EQ(3u, sp[4]);
    ASSERT_EQ(3u, sp[5]);
    ASSERT_EQ(3u, sp[6]);
}


TEST(NumLib_SparsityPattern, MultipleComponentsLinearMesh)
{
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::MeshGenerator::generateLineMesh(3u, 1.));
    std::unique_ptr<MeshLib::MeshSubset const> nodesSubset(
        new MeshLib::MeshSubset(*mesh, &mesh->getNodes()));

    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> components;
    components.emplace_back(new MeshLib::MeshSubsets{nodesSubset.get()});
    components.emplace_back(new MeshLib::MeshSubsets{nodesSubset.get()});
    NumLib::LocalToGlobalIndexMap dof_map(
                      std::move(components),
                      NumLib::ComponentOrder::BY_COMPONENT);

    GlobalSparsityPattern sp = NumLib::computeSparsityPattern(dof_map, *mesh.get());

    ASSERT_EQ(8u, sp.size());
    for (int i=0; i<2; i++)
    {
        ASSERT_EQ(4u, sp[i*mesh->getNumberOfNodes() + 0]);
        ASSERT_EQ(6u, sp[i*mesh->getNumberOfNodes() + 1]);
        ASSERT_EQ(6u, sp[i*mesh->getNumberOfNodes() + 2]);
        ASSERT_EQ(4u, sp[i*mesh->getNumberOfNodes() + 3]);
    }
}


TEST(NumLib_SparsityPattern, MultipleComponentsLinearQuadraticMesh)
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

    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> components;
    components.emplace_back(new MeshLib::MeshSubsets{baseNodesSubset.get()});
    components.emplace_back(new MeshLib::MeshSubsets{allNodesSubset.get()});
    NumLib::LocalToGlobalIndexMap dof_map(
                      std::move(components),
                      NumLib::ComponentOrder::BY_COMPONENT);

    GlobalSparsityPattern sp = NumLib::computeSparsityPattern(dof_map, *mesh.get());

    ASSERT_EQ(11u, sp.size());
    // 1st component
    ASSERT_EQ(5u, sp[0]);
    ASSERT_EQ(8u, sp[1]);
    ASSERT_EQ(8u, sp[2]);
    ASSERT_EQ(5u, sp[3]);
    // 2nd component
    ASSERT_EQ(5u, sp[4]);
    ASSERT_EQ(8u, sp[5]);
    ASSERT_EQ(8u, sp[6]);
    ASSERT_EQ(5u, sp[7]);
    ASSERT_EQ(5u, sp[8]);
    ASSERT_EQ(5u, sp[9]);
    ASSERT_EQ(5u, sp[10]);
}

