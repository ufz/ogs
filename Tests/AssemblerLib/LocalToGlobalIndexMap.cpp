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

#include "MeshLib/MeshSearcher.h"
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

    // There must be as many rows as nodes in the input times the number of
    // components.
    ASSERT_EQ(mesh->getNNodes() * components.size(), dof_map->dofSize());
}

TEST_F(AssemblerLibLocalToGlobalIndexMapTest, NumberOfRowsByLocation)
{
    dof_map = new AssemblerLib::LocalToGlobalIndexMap(components,
        AssemblerLib::ComponentOrder::BY_LOCATION);

    // There must be as many rows as nodes in the input times the number of
    // components.
    ASSERT_EQ(mesh->getNNodes() * components.size(), dof_map->dofSize());
}

TEST_F(AssemblerLibLocalToGlobalIndexMapTest, SubsetByComponent)
{
    dof_map = new AssemblerLib::LocalToGlobalIndexMap(components,
        AssemblerLib::ComponentOrder::BY_COMPONENT);

    // Select some elements from the full mesh.
    std::array<std::size_t, 3> const ids = {{ 0, 5, 8 }};
    std::vector<MeshLib::Element*> some_elements;
    for (std::size_t id : ids)
        some_elements.push_back(const_cast<MeshLib::Element*>(mesh->getElement(id)));

    // Find unique node ids of the selected elements for testing.
    std::vector<MeshLib::Node*> selected_nodes = MeshLib::selectNodes(some_elements);

    MeshLib::MeshSubset const* const selected_subset =
        nodesSubset->getIntersectionByNodes(selected_nodes);
    std::vector<MeshLib::MeshSubsets*> selected_components;
    selected_components.emplace_back(new MeshLib::MeshSubsets(selected_subset));
    selected_components.emplace_back(new MeshLib::MeshSubsets(selected_subset));

    AssemblerLib::LocalToGlobalIndexMap* dof_map_subset =
        dof_map->deriveBoundaryConstrainedMap(selected_components, some_elements);

    // There must be as many rows as nodes in the input times the number of
    // components.
    ASSERT_EQ(selected_nodes.size() * selected_components.size(),
            dof_map_subset->dofSize());
}
