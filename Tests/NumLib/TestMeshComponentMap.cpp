/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "NumLib/DOF/MeshComponentMap.h"

#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshSubsets.h"

class NumLibMeshComponentMapTest : public ::testing::Test
{
    public:
    typedef MeshLib::MeshItemType MeshItemType;
    typedef MeshLib::Location Location;
    typedef NumLib::MeshComponentMap MeshComponentMap;

    public:
    NumLibMeshComponentMapTest()
        : mesh(nullptr), nodesSubset(nullptr), cmap(nullptr)
    {
        mesh = MeshLib::MeshGenerator::generateLineMesh(1.0, mesh_size);
        nodesSubset = new MeshLib::MeshSubset(*mesh, &mesh->getNodes());

        // Add two components both based on the same nodesSubset.
        components.emplace_back(new MeshLib::MeshSubsets{nodesSubset});
        components.emplace_back(new MeshLib::MeshSubsets{nodesSubset});
    }

    ~NumLibMeshComponentMapTest()
    {
        delete cmap;
        delete nodesSubset;
        delete mesh;
    }

    static std::size_t const mesh_size = 9;
    MeshLib::Mesh const* mesh;
    MeshLib::MeshSubset const* nodesSubset;

    //data component 0 and 1 are assigned to all nodes in the mesh
    static std::size_t const comp0_id = 0;
    static std::size_t const comp1_id = 1;
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> components;
    MeshComponentMap const* cmap;

    //
    // Functions used for checking.
    //

    // Returns global index of a node location and a component.
    std::size_t giAtNodeForComponent(std::size_t const n, std::size_t const c) const
    {
        return cmap->getGlobalIndex(Location(mesh->getID(), MeshItemType::Node, n), c);
    }

};

#ifndef USE_PETSC
TEST_F(NumLibMeshComponentMapTest, CheckOrderByComponent)
#else
TEST_F(NumLibMeshComponentMapTest, DISABLED_CheckOrderByComponent)
#endif
{
    // - Entries in the vector are arranged in the order of a component type and then node ID
    // - For example, x=[(node 0, comp 0) (node 1, comp 0) ... (node n, comp0), (node 0, comp1) ... ]

    cmap = new MeshComponentMap(components,
        NumLib::ComponentOrder::BY_COMPONENT);

    ASSERT_EQ(2 * mesh->getNumberOfNodes(), cmap->dofSizeWithGhosts());
    for (std::size_t i = 0; i < mesh_size; i++)
    {
        // Test global indices for the different components of the node.
        ASSERT_EQ(i , giAtNodeForComponent(i, comp0_id));
        ASSERT_EQ(mesh_size + 1 + i, giAtNodeForComponent(i, comp1_id));

        // Test component ids of the node.
        std::vector<std::size_t> const vecCompIDs = cmap->getComponentIDs(
            Location(mesh->getID(), MeshItemType::Node, i));
        ASSERT_EQ(2u, vecCompIDs.size());
        ASSERT_EQ(0u, vecCompIDs[0]);
        ASSERT_EQ(1u, vecCompIDs[1]);
    }
}

#ifndef USE_PETSC
TEST_F(NumLibMeshComponentMapTest, CheckOrderByLocation)
#else
TEST_F(NumLibMeshComponentMapTest, DISABLED_CheckOrderByLocation)
#endif
{
    // - Entries in the vector are arranged in the order of node ID and then a component type
    // - For example, x=[(node 0, comp 0) (node 0, comp 1) ... (node n, comp0), (node n, comp1) ]

    cmap = new MeshComponentMap(components,
        NumLib::ComponentOrder::BY_LOCATION);

    ASSERT_EQ(2 * mesh->getNumberOfNodes(), cmap->dofSizeWithGhosts());
    for (std::size_t i = 0; i < mesh_size; i++)
    {
        // Test global indices for the different components of the node.
        ASSERT_EQ(2 * i , giAtNodeForComponent(i, comp0_id));
        ASSERT_EQ(2 * i + 1, giAtNodeForComponent(i, comp1_id));

        // Test component ids of the node.
        std::vector<std::size_t> const vecCompIDs = cmap->getComponentIDs(
            Location(mesh->getID(), MeshItemType::Node, i));
        ASSERT_EQ(2u, vecCompIDs.size());
        ASSERT_EQ(0u, vecCompIDs[0]);
        ASSERT_EQ(1u, vecCompIDs[1]);
    }
}

#ifndef USE_PETSC
TEST_F(NumLibMeshComponentMapTest, OutOfRangeAccess)
#else
TEST_F(NumLibMeshComponentMapTest, DISABLED_OutOfRangeAccess)
#endif
{
    cmap = new MeshComponentMap(components,
        NumLib::ComponentOrder::BY_COMPONENT);

    ASSERT_EQ(MeshComponentMap::nop, cmap->getGlobalIndex(
        Location(mesh->getID(), MeshItemType::Node, mesh_size + 1), comp0_id));
    ASSERT_EQ(MeshComponentMap::nop, cmap->getGlobalIndex(
        Location(mesh->getID() + 1, MeshItemType::Node, 0), comp0_id));
    ASSERT_EQ(MeshComponentMap::nop, cmap->getGlobalIndex(
        Location(mesh->getID(), MeshItemType::Cell, 0), comp0_id));
    ASSERT_EQ(MeshComponentMap::nop, cmap->getGlobalIndex(
        Location(mesh->getID(), MeshItemType::Node, 0), 10));
}

#ifndef USE_PETSC
TEST_F(NumLibMeshComponentMapTest, SubsetOfNodesByComponent)
#else
TEST_F(NumLibMeshComponentMapTest, DISABLED_SubsetOfNodesByComponent)
#endif
{
    cmap = new MeshComponentMap(components,
        NumLib::ComponentOrder::BY_COMPONENT);

    // Select some nodes from the full mesh.
    std::array<std::size_t, 3> const ids = {{ 0, 5, 9 }};
    std::vector<MeshLib::Node*> some_nodes;
    for (std::size_t id : ids)
        some_nodes.push_back(const_cast<MeshLib::Node*>(mesh->getNode(id)));

    MeshLib::MeshSubset some_nodes_mesh_subset(*mesh, &some_nodes);

    std::size_t const selected_component_id = 1;
    auto selected_component = std::unique_ptr<MeshLib::MeshSubsets>{
        new MeshLib::MeshSubsets{&some_nodes_mesh_subset}};

    // Subset the original cmap.
    MeshComponentMap cmap_subset =
        cmap->getSubset(selected_component_id, *selected_component);

    // Check number of components as selected
    ASSERT_EQ(ids.size(), cmap_subset.dofSizeWithGhosts());

    // .. and the content of the subset.
    for (std::size_t id : ids)
    {
        Location const l(mesh->getID(), MeshItemType::Node, id);
        EXPECT_EQ(cmap->getGlobalIndex(l, comp1_id),
            cmap_subset.getGlobalIndex(l, comp1_id));
    }
}

#ifndef USE_PETSC
TEST_F(NumLibMeshComponentMapTest, SubsetOfNodesByLocation)
#else
TEST_F(NumLibMeshComponentMapTest, DISABLED_SubsetOfNodesByLocation)
#endif
{
    cmap = new MeshComponentMap(components,
        NumLib::ComponentOrder::BY_LOCATION);

    // Select some nodes from the full mesh.
    std::array<std::size_t, 3> const ids = {{ 0, 5, 9 }};
    std::vector<MeshLib::Node*> some_nodes;
    for (std::size_t id : ids)
        some_nodes.push_back(const_cast<MeshLib::Node*>(mesh->getNode(id)));

    MeshLib::MeshSubset some_nodes_mesh_subset(*mesh, &some_nodes);

    std::size_t const selected_component_id = 1;
    auto selected_component = std::unique_ptr<MeshLib::MeshSubsets>{
        new MeshLib::MeshSubsets{&some_nodes_mesh_subset}};

    // Subset the original cmap.
    MeshComponentMap cmap_subset =
        cmap->getSubset(selected_component_id, *selected_component);

    // Check number of components as selected
    ASSERT_EQ(ids.size(), cmap_subset.dofSizeWithGhosts());

    // .. and the content of the subset.
    for (std::size_t id : ids)
    {
        Location const l(mesh->getID(), MeshItemType::Node, id);
        EXPECT_EQ(cmap->getGlobalIndex(l, comp1_id),
            cmap_subset.getGlobalIndex(l, comp1_id));
    }
}
