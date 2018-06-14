/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "MeshLib/MeshSubset.h"

class NumLibMeshComponentMapTest : public ::testing::Test
{
    public:
        using MeshItemType = MeshLib::MeshItemType;
        using Location = MeshLib::Location;
        using MeshComponentMap = NumLib::MeshComponentMap;

    public:
        NumLibMeshComponentMapTest() : mesh(nullptr), cmap(nullptr)
        {
            mesh = MeshLib::MeshGenerator::generateLineMesh(1.0, mesh_size);
            MeshLib::MeshSubset nodesSubset{*mesh, mesh->getNodes()};

            // Add two components both based on the same nodesSubset.
            components.emplace_back(nodesSubset);
            components.emplace_back(nodesSubset);
    }

    ~NumLibMeshComponentMapTest() override
    {
        delete cmap;
        delete mesh;
    }

    static std::size_t const mesh_size = 9;
    MeshLib::Mesh const* mesh;

    //data component 0 and 1 are assigned to all nodes in the mesh
    static std::size_t const comp0_id = 0;
    static std::size_t const comp1_id = 1;
    std::vector<MeshLib::MeshSubset> components;
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
        std::vector<int> const vecCompIDs = cmap->getComponentIDs(
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
        std::vector<int> const vecCompIDs = cmap->getComponentIDs(
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

MeshLib::Mesh createMeshFromSelectedNodes(
    MeshLib::Mesh const& mesh, std::vector<std::size_t> const& selected_nodes)
{
    // Deep copy of selected nodes from the mesh.
    std::vector<MeshLib::Node*> some_nodes;
    std::transform(begin(selected_nodes), end(selected_nodes),
                   back_inserter(some_nodes),
                   [&mesh](std::size_t const node_id) {
                       return new MeshLib::Node(*mesh.getNode(node_id));
                   });

    // The resulting mesh without elements containing the selected nodes.
    MeshLib::Mesh result("boundary_mesh", some_nodes, {});
    addPropertyToMesh(result, "bulk_node_ids", MeshLib::MeshItemType::Node, 1,
                      selected_nodes);
    return result;
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
    std::vector<std::size_t> const ids = {0, 5, 9};
    // A smaller mesh without elements containing the selected nodes.
    auto boundary_mesh = createMeshFromSelectedNodes(*mesh, ids);

    MeshLib::MeshSubset const selected_component(boundary_mesh,
                                                 boundary_mesh.getNodes());

    int const selected_component_id = 1;

    // Subset the original cmap.
    MeshComponentMap const cmap_subset = cmap->getSubset(
        components, selected_component, {selected_component_id});

    // Check number of components as selected
    ASSERT_EQ(ids.size(), cmap_subset.dofSizeWithGhosts());

    // .. and the content of the subset.
    for (auto const* n : boundary_mesh.getNodes())
    {
        std::size_t const id = n->getID();
        Location const l_bulk(mesh->getID(), MeshItemType::Node, ids[id]);
        Location const l_boundary(boundary_mesh.getID(), MeshItemType::Node,
                                  id);
        EXPECT_EQ(cmap->getGlobalIndex(l_bulk, comp1_id),
                  cmap_subset.getGlobalIndex(l_boundary, comp1_id));
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
    std::vector<std::size_t> const ids = {0, 5, 9};
    // A smaller mesh without elements containing the selected nodes.
    auto boundary_mesh = createMeshFromSelectedNodes(*mesh, ids);

    MeshLib::MeshSubset const selected_component(boundary_mesh,
                                                 boundary_mesh.getNodes());

    int const selected_component_id = 1;

    // Subset the original cmap.
    MeshComponentMap const cmap_subset = cmap->getSubset(
        components, selected_component, {selected_component_id});

    // Check number of components as selected
    ASSERT_EQ(ids.size(), cmap_subset.dofSizeWithGhosts());

    // .. and the content of the subset.
    for (auto const* n : boundary_mesh.getNodes())
    {
        std::size_t const id = n->getID();
        Location const l_bulk(mesh->getID(), MeshItemType::Node, ids[id]);
        Location const l_boundary(boundary_mesh.getID(), MeshItemType::Node,
                                  id);
        EXPECT_EQ(cmap->getGlobalIndex(l_bulk, comp1_id),
                  cmap_subset.getGlobalIndex(l_boundary, comp1_id));
    }
}

#ifndef USE_PETSC
TEST_F(NumLibMeshComponentMapTest, MulticomponentVariable)
#else
TEST_F(NumLibMeshComponentMapTest, DISABLED_MulticomponentVariable)
#endif
{
    cmap =
        new MeshComponentMap(components, NumLib::ComponentOrder::BY_LOCATION);

    // Select some nodes from the full mesh.
    std::vector<std::size_t> const ids = {0, 5, 9};
    // A smaller mesh without elements containing the selected nodes.
    auto boundary_mesh = createMeshFromSelectedNodes(*mesh, ids);

    MeshLib::MeshSubset const selected_component(boundary_mesh,
                                                 boundary_mesh.getNodes());

    // Subset the original cmap.
    std::vector<int> const selected_component_ids = {0, 1};
    MeshComponentMap const cmap_subset =
        cmap->getSubset(components, selected_component, selected_component_ids);

    // Check number of components as selected
    ASSERT_EQ(ids.size() * selected_component_ids.size(),
              cmap_subset.dofSizeWithGhosts());

    // .. and the content of the subset.
    for (auto const* n : boundary_mesh.getNodes())
    {
        std::size_t const id = n->getID();
        Location const l_bulk(mesh->getID(), MeshItemType::Node, ids[id]);
        Location const l_boundary(boundary_mesh.getID(), MeshItemType::Node,
                                  id);
        for (auto const& c : selected_component_ids)
            EXPECT_EQ(cmap->getGlobalIndex(l_bulk, c),
                      cmap_subset.getGlobalIndex(l_boundary, c));
    }
}

#ifndef USE_PETSC
TEST_F(NumLibMeshComponentMapTest, MulticomponentVariableSingleComponent)
#else
TEST_F(NumLibMeshComponentMapTest,
       DISABLED_MulticomponentVariableSingleComponent)
#endif
{
    cmap =
        new MeshComponentMap(components, NumLib::ComponentOrder::BY_LOCATION);

    // Select some nodes from the full mesh.
    std::vector<std::size_t> const ids = {0, 5, 9};
    // A smaller mesh without elements containing the selected nodes.
    auto boundary_mesh = createMeshFromSelectedNodes(*mesh, ids);

    MeshLib::MeshSubset const selected_component(boundary_mesh,
                                                 boundary_mesh.getNodes());

    // Subset the original cmap.
    std::vector<int> const selected_component_ids = {1};
    MeshComponentMap const cmap_subset =
        cmap->getSubset(components, selected_component, selected_component_ids);

    // Check number of components as selected
    ASSERT_EQ(ids.size() * selected_component_ids.size(),
              cmap_subset.dofSizeWithGhosts());

    // .. and the content of the subset.
    for (auto const* n : boundary_mesh.getNodes())
    {
        std::size_t const id = n->getID();
        Location const l_bulk(mesh->getID(), MeshItemType::Node, ids[id]);
        Location const l_boundary(boundary_mesh.getID(), MeshItemType::Node,
                                  id);
        for (auto const& c : selected_component_ids)
            EXPECT_EQ(cmap->getGlobalIndex(l_bulk, c),
                      cmap_subset.getGlobalIndex(l_boundary, c));
    }
}
