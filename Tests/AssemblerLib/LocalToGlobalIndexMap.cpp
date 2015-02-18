/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <vector>

#include <gtest/gtest.h>

#include "AssemblerLib/LocalToGlobalIndexMap.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"

#include "MeshLib/MeshSubsets.h"
#include "MeshLib/Mesh.h"


class AssemblerLibLocalToGlobalIndexMapTest : public ::testing::Test
{

public:
    AssemblerLibLocalToGlobalIndexMapTest()
    {
        mesh = MeshLib::MeshGenerator::generateLineMesh(1.0, mesh_size);
        nodesSubset = new MeshLib::MeshSubset(*mesh, mesh->getNodes());

        // Add two components both based on the same nodesSubset.
        components.emplace_back(new MeshLib::MeshSubsets(nodesSubset));
        components.emplace_back(new MeshLib::MeshSubsets(nodesSubset));
    }

    ~AssemblerLibLocalToGlobalIndexMapTest()
    {
        delete dof_map;
        for (auto p : components)
            delete p;

        delete nodesSubset;
        delete mesh;
    }

protected:
    static std::size_t const mesh_size = 9;
    MeshLib::Mesh const* mesh = nullptr;
    MeshLib::MeshSubset const* nodesSubset = nullptr;

    //data component 0 and 1 are assigned to all nodes in the mesh
    static std::size_t const comp0_id = 0;
    static std::size_t const comp1_id = 1;
    std::vector<MeshLib::MeshSubsets*> components;

    AssemblerLib::LocalToGlobalIndexMap const* dof_map = nullptr;
};

TEST_F(AssemblerLibLocalToGlobalIndexMapTest, NumberOfRowsByComponent)
{
    dof_map = new AssemblerLib::LocalToGlobalIndexMap(components,
        AssemblerLib::ComponentOrder::BY_COMPONENT);

    //std::cout << "# dof table\n" << *dof_map << std::endl;

    // There must be as many rows as nodes in the input times the number of
    // components.
    ASSERT_EQ(mesh->getNNodes() * components.size(), dof_map->dofSize());
}

TEST_F(AssemblerLibLocalToGlobalIndexMapTest, NumberOfRowsByLocation)
{
    dof_map = new AssemblerLib::LocalToGlobalIndexMap(components,
        AssemblerLib::ComponentOrder::BY_LOCATION);

    //std::cout << "# dof table\n" << *dof_map << std::endl;

    // There must be as many rows as nodes in the input times the number of
    // components.
    ASSERT_EQ(mesh->getNNodes() * components.size(), dof_map->dofSize());
}
