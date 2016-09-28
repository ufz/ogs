/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <memory>
#include <vector>

#include <gtest/gtest.h>

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"

#include "MeshLib/MeshSearch/NodeSearch.h"
#include "MeshLib/MeshSubsets.h"
#include "MeshLib/Mesh.h"


class NumLibLocalToGlobalIndexMapTest : public ::testing::Test
{

public:
    NumLibLocalToGlobalIndexMapTest()
    {
        mesh.reset(MeshLib::MeshGenerator::generateLineMesh(1.0, mesh_size));
        nodesSubset.reset(new MeshLib::MeshSubset(*mesh, &mesh->getNodes()));

        // Add two components both based on the same nodesSubset.
        components.emplace_back(new MeshLib::MeshSubsets{nodesSubset.get()});
        components.emplace_back(new MeshLib::MeshSubsets{nodesSubset.get()});
    }

protected:
    static std::size_t const mesh_size = 9;
    std::unique_ptr<MeshLib::Mesh const> mesh;
    std::unique_ptr<MeshLib::MeshSubset const> nodesSubset;

    //data component 0 and 1 are assigned to all nodes in the mesh
    static std::size_t const comp0_id = 0;
    static std::size_t const comp1_id = 1;
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> components;

    std::unique_ptr<NumLib::LocalToGlobalIndexMap const> dof_map;
};

#ifndef USE_PETSC
TEST_F(NumLibLocalToGlobalIndexMapTest, NumberOfRowsByComponent)
#else
TEST_F(NumLibLocalToGlobalIndexMapTest, DISABLED_NumberOfRowsByComponent)
#endif
{
    // need to store the size because the components will be moved into the
    // DOF-table.
    std::size_t components_size = components.size();

    dof_map.reset(new NumLib::LocalToGlobalIndexMap(std::move(components),
        NumLib::ComponentOrder::BY_COMPONENT));

    // There must be as many rows as nodes in the input times the number of
    // components.
    ASSERT_EQ(mesh->getNumberOfNodes() * components_size, dof_map->dofSizeWithGhosts());
}

#ifndef USE_PETSC
TEST_F(NumLibLocalToGlobalIndexMapTest, NumberOfRowsByLocation)
#else
TEST_F(NumLibLocalToGlobalIndexMapTest, DISABLED_NumberOfRowsByLocation)
#endif
{
    // need to store the size because the components will be moved into the
    // DOF-table.
    std::size_t components_size = components.size();

    dof_map.reset(new NumLib::LocalToGlobalIndexMap(std::move(components),
        NumLib::ComponentOrder::BY_LOCATION));

    // There must be as many rows as nodes in the input times the number of
    // components.
    ASSERT_EQ(mesh->getNumberOfNodes() * components_size, dof_map->dofSizeWithGhosts());
}

#ifndef USE_PETSC
TEST_F(NumLibLocalToGlobalIndexMapTest, SubsetByComponent)
#else
TEST_F(NumLibLocalToGlobalIndexMapTest, DISABLED_SubsetByComponent)
#endif
{
    dof_map.reset(new NumLib::LocalToGlobalIndexMap(std::move(components),
        NumLib::ComponentOrder::BY_COMPONENT));

    // Select some elements from the full mesh.
    std::array<std::size_t, 3> const ids = {{ 0, 5, 8 }};
    std::vector<MeshLib::Element*> some_elements;
    for (std::size_t id : ids)
        some_elements.push_back(const_cast<MeshLib::Element*>(mesh->getElement(id)));

    // Find unique node ids of the selected elements for testing.
    std::vector<MeshLib::Node*> selected_nodes = MeshLib::getUniqueNodes(some_elements);

    auto const selected_subset = std::unique_ptr<MeshLib::MeshSubset const>{
        nodesSubset->getIntersectionByNodes(selected_nodes)};
    auto selected_component = std::unique_ptr<MeshLib::MeshSubsets>{
        new MeshLib::MeshSubsets{selected_subset.get()}};

    auto dof_map_subset = std::unique_ptr<NumLib::LocalToGlobalIndexMap>{
        dof_map->deriveBoundaryConstrainedMap(1,  // variable id
                                              0,  // component id
                                              std::move(selected_component),
                                              some_elements)};

    // There must be as many rows as nodes in the input times the number of
    // components.
    ASSERT_EQ(selected_nodes.size(), dof_map_subset->dofSizeWithGhosts());
}
