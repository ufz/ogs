/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>

#include "gtest/gtest.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSubsets.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"

using namespace MeshLib;

TEST(MeshLibMeshSubsets, UniqueMeshIds)
{
    // Create first mesh
    Mesh const m0("first", std::vector<Node*>(), std::vector<Element*>());
    Mesh const m1("second", std::vector<Node*>(), std::vector<Element*>());

    std::vector<Node*> const empty_node_ptr_vector(0);

    MeshSubset const ms0(m0, &empty_node_ptr_vector);
    MeshSubset const ms1(m1, &empty_node_ptr_vector);
    MeshSubset const ms1a(m1, &empty_node_ptr_vector);

    MeshSubset const* const mesh_subsets[3] = {&ms0, &ms1, &ms1a};

    // EXPECT_NO_DEATH
    MeshSubsets(&mesh_subsets[0], &mesh_subsets[0] + 2);

    EXPECT_THROW(MeshSubsets(&mesh_subsets[1], &mesh_subsets[1] + 2), std::runtime_error);
    EXPECT_THROW(MeshSubsets(&mesh_subsets[0], &mesh_subsets[0] + 3), std::runtime_error);
}


TEST(MeshLibMeshSubsets, GetIntersectionByNodes)
{
    auto mesh = std::unique_ptr<Mesh>{MeshGenerator::generateLineMesh(1., 10)};
    MeshSubset all_nodes_mesh_subset(*mesh, &mesh->getNodes());

    // Select nodes
    std::vector<Node*> some_nodes;
    some_nodes.assign({
        const_cast<Node*>(mesh->getNode(0)),
        const_cast<Node*>(mesh->getNode(2)),
        const_cast<Node*>(mesh->getNode(5)),
        const_cast<Node*>(mesh->getNode(7)) });
    auto some_nodes_mesh_subset = std::unique_ptr<MeshSubset>
        {all_nodes_mesh_subset.getIntersectionByNodes(some_nodes)};

    // Check sizes.
    ASSERT_EQ(some_nodes.size(), some_nodes_mesh_subset->getNumberOfNodes());

    // Check ids.
    std::size_t nnodes = some_nodes.size();
    for (std::size_t i = 0; i < nnodes; ++i)
        ASSERT_EQ(some_nodes[i]->getID(), some_nodes_mesh_subset->getNodeID(i));
}

