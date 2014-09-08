/**
 * \author Norihiro Watanabe
 * \date   2013-04-16
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <vector>

#include "AssemblerLib/MeshComponentMap.h"

#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshSubsets.h"

class AssemblerLibMeshComponentMapTest : public ::testing::Test
{
    public:
    typedef MeshLib::MeshItemType MeshItemType;
    typedef MeshLib::Location Location;
    typedef AssemblerLib::MeshComponentMap MeshComponentMap;

    public:
    AssemblerLibMeshComponentMapTest()
        : mesh(nullptr), nodesSubset(nullptr), cmap(nullptr)
    {
        mesh = MeshLib::MeshGenerator::generateLineMesh(1.0, mesh_size);
        nodesSubset = new MeshLib::MeshSubset(*mesh, mesh->getNodes());

        // Add two components both based on the same nodesSubset.
        components.emplace_back(new MeshLib::MeshSubsets(nodesSubset));
        components.emplace_back(new MeshLib::MeshSubsets(nodesSubset));
    }

    ~AssemblerLibMeshComponentMapTest()
    {
        delete cmap;
        std::remove_if(components.begin(), components.end(),
            [](MeshLib::MeshSubsets* p) { delete p; return true; });
        delete nodesSubset;
        delete mesh;
    }

    static std::size_t const mesh_size = 9;
    MeshLib::Mesh const* mesh;
    MeshLib::MeshSubset const* nodesSubset;

    //data component 0 and 1 are assigned to all nodes in the mesh
    static std::size_t const comp0_id = 0;
    static std::size_t const comp1_id = 1;
    std::vector<MeshLib::MeshSubsets*> components;
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

TEST_F(AssemblerLibMeshComponentMapTest, CheckOrderByComponent)
{
    // - Entries in the vector are arranged in the order of a component type and then node ID
    // - For example, x=[(node 0, comp 0) (node 1, comp 0) ... (node n, comp0), (node 0, comp1) ... ]

    cmap = new MeshComponentMap(components,
        AssemblerLib::ComponentOrder::BY_COMPONENT);

    //std::cout << "# database \n" << *cmap << std::endl;

    ASSERT_EQ(2 * mesh->getNNodes(), cmap->size());
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

TEST_F(AssemblerLibMeshComponentMapTest, CheckOrderByLocation)
{
    // - Entries in the vector are arranged in the order of node ID and then a component type
    // - For example, x=[(node 0, comp 0) (node 0, comp 1) ... (node n, comp0), (node n, comp1) ]

    cmap = new MeshComponentMap(components,
        AssemblerLib::ComponentOrder::BY_LOCATION);

    //std::cout << "# database \n" << *cmap << std::endl;

    ASSERT_EQ(2 * mesh->getNNodes(), cmap->size());
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

TEST_F(AssemblerLibMeshComponentMapTest, OutOfRangeAccess)
{
    cmap = new MeshComponentMap(components,
        AssemblerLib::ComponentOrder::BY_COMPONENT);

    //std::cout << "# database \n" << *cmap << std::endl;

    ASSERT_EQ(MeshComponentMap::nop, cmap->getGlobalIndex(
        Location(mesh->getID(), MeshItemType::Node, mesh_size + 1), comp0_id));
    ASSERT_EQ(MeshComponentMap::nop, cmap->getGlobalIndex(
        Location(mesh->getID() + 1, MeshItemType::Node, 0), comp0_id));
    ASSERT_EQ(MeshComponentMap::nop, cmap->getGlobalIndex(
        Location(mesh->getID(), MeshItemType::Cell, 0), comp0_id));
    ASSERT_EQ(MeshComponentMap::nop, cmap->getGlobalIndex(
        Location(mesh->getID(), MeshItemType::Node, 0), 10));
}
