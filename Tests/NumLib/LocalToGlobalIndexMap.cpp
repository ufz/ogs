/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
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
        nodesSubset =
            std::make_unique<MeshLib::MeshSubset>(*mesh, &mesh->getNodes());

        // Add two components both based on the same nodesSubset.
        components.emplace_back(nodesSubset.get());
        components.emplace_back(nodesSubset.get());
    }

protected:
    static std::size_t const mesh_size = 9;
    std::unique_ptr<MeshLib::Mesh const> mesh;
    std::unique_ptr<MeshLib::MeshSubset const> nodesSubset;

    //data component 0 and 1 are assigned to all nodes in the mesh
    static int const comp0_id = 0;
    static int const comp1_id = 1;
    std::vector<MeshLib::MeshSubsets> components;

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

    dof_map = std::make_unique<NumLib::LocalToGlobalIndexMap>(
        std::move(components), NumLib::ComponentOrder::BY_COMPONENT);

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

    dof_map = std::make_unique<NumLib::LocalToGlobalIndexMap>(
        std::move(components), NumLib::ComponentOrder::BY_LOCATION);

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
    dof_map = std::make_unique<NumLib::LocalToGlobalIndexMap>(
        std::move(components), NumLib::ComponentOrder::BY_COMPONENT);

    // Select some elements from the full mesh.
    std::array<std::size_t, 3> const ids = {{ 0, 5, 8 }};
    std::vector<MeshLib::Element*> some_elements;
    for (std::size_t id : ids)
        some_elements.push_back(const_cast<MeshLib::Element*>(mesh->getElement(id)));

    // Find unique node ids of the selected elements for testing.
    std::vector<MeshLib::Node*> selected_nodes = MeshLib::getUniqueNodes(some_elements);

    auto const selected_subset = std::unique_ptr<MeshLib::MeshSubset const>{
        nodesSubset->getIntersectionByNodes(selected_nodes)};
    MeshLib::MeshSubsets selected_component{selected_subset.get()};

    auto dof_map_subset = std::unique_ptr<NumLib::LocalToGlobalIndexMap>{
        dof_map->deriveBoundaryConstrainedMap(1,    // variable id
                                              {0},  // component id
                                              std::move(selected_component),
                                              some_elements)};

    // There must be as many rows as nodes in the input times the number of
    // components.
    ASSERT_EQ(selected_nodes.size(), dof_map_subset->dofSizeWithGhosts());
}

#ifndef USE_PETSC
TEST_F(NumLibLocalToGlobalIndexMapTest, MultipleVariablesMultipleComponents)
#else
TEST_F(NumLibLocalToGlobalIndexMapTest, DISABLED_MultipleVariablesMultipleComponents)
#endif
{
    // test 2 variables (1st variable with 1 component, 2nd variable with 2 components)
    components.emplace_back(nodesSubset.get());

    std::vector<int> vec_var_n_components{1, 2};

    dof_map = std::make_unique<NumLib::LocalToGlobalIndexMap>(
        std::move(components),
        vec_var_n_components,
        NumLib::ComponentOrder::BY_COMPONENT);

    ASSERT_EQ(30, dof_map->dofSizeWithGhosts());
    ASSERT_EQ(3, dof_map->getNumberOfComponents());
    ASSERT_EQ(2u, dof_map->getNumberOfVariables());
    ASSERT_EQ(1, dof_map->getNumberOfVariableComponents(0));
    ASSERT_EQ(2, dof_map->getNumberOfVariableComponents(1));

    MeshLib::Location l_node0(mesh->getID(), MeshLib::MeshItemType::Node, 0);

    ASSERT_EQ(0, dof_map->getGlobalIndex(l_node0, 0, 0));
    ASSERT_EQ(10, dof_map->getGlobalIndex(l_node0, 1, 0));
    ASSERT_EQ(20, dof_map->getGlobalIndex(l_node0, 1, 1));
}

#ifndef USE_PETSC
TEST_F(NumLibLocalToGlobalIndexMapTest, MultipleVariablesMultipleComponents2)
#else
TEST_F(NumLibLocalToGlobalIndexMapTest, DISABLED_MultipleVariablesMultipleComponents2)
#endif
{
    // test 2 variables (1st variable with 2 component, 2nd variable with 1 components)
    components.emplace_back(nodesSubset.get());

    std::vector<int> vec_var_n_components{2, 1};

    dof_map = std::make_unique<NumLib::LocalToGlobalIndexMap>(
        std::move(components),
        vec_var_n_components,
        NumLib::ComponentOrder::BY_COMPONENT);

    ASSERT_EQ(30, dof_map->dofSizeWithGhosts());
    ASSERT_EQ(3, dof_map->getNumberOfComponents());
    ASSERT_EQ(2u, dof_map->getNumberOfVariables());
    ASSERT_EQ(2, dof_map->getNumberOfVariableComponents(0));
    ASSERT_EQ(1, dof_map->getNumberOfVariableComponents(1));

    MeshLib::Location l_node0(mesh->getID(), MeshLib::MeshItemType::Node, 0);

    ASSERT_EQ(0, dof_map->getGlobalIndex(l_node0, 0, 0));
    ASSERT_EQ(10, dof_map->getGlobalIndex(l_node0, 0, 1));
    ASSERT_EQ(20, dof_map->getGlobalIndex(l_node0, 1, 0));
}


#ifndef USE_PETSC
TEST_F(NumLibLocalToGlobalIndexMapTest, MultipleVariablesMultipleComponentsHeterogeneousElements)
#else
TEST_F(NumLibLocalToGlobalIndexMapTest, DISABLED_MultipleVariablesMultipleComponentsHeterogeneousElements)
#endif
{
    // test 2 variables
    // - 1st variable with 2 components for all nodes, elements
    // - 2nd variable with 1 component for nodes of element id 1
    std::vector<MeshLib::Node*> var2_nodes{const_cast<MeshLib::Node*>(mesh->getNode(1)), const_cast<MeshLib::Node*>(mesh->getNode(2))};
    auto var2_subset =
        std::make_unique<MeshLib::MeshSubset>(*mesh, &var2_nodes);
    components.emplace_back(var2_subset.get());

    std::vector<int> vec_var_n_components{2, 1};
    std::vector<std::vector<MeshLib::Element*>const*> vec_var_elements;
    vec_var_elements.push_back(&mesh->getElements());
    std::vector<MeshLib::Element*> var2_elements{const_cast<MeshLib::Element*>(mesh->getElement(1))};
    vec_var_elements.push_back(&var2_elements);

    dof_map = std::make_unique<NumLib::LocalToGlobalIndexMap>(
        std::move(components),
        vec_var_n_components,
        vec_var_elements,
        NumLib::ComponentOrder::BY_COMPONENT);

    ASSERT_EQ(22u, dof_map->dofSizeWithGhosts());
    ASSERT_EQ(3, dof_map->getNumberOfComponents());
    ASSERT_EQ(2u, dof_map->getNumberOfVariables());
    ASSERT_EQ(2, dof_map->getNumberOfVariableComponents(0));
    ASSERT_EQ(1, dof_map->getNumberOfVariableComponents(1));

    MeshLib::Location l_node0(mesh->getID(), MeshLib::MeshItemType::Node, 0);

    ASSERT_EQ(0, dof_map->getGlobalIndex(l_node0, 0, 0));
    ASSERT_EQ(10, dof_map->getGlobalIndex(l_node0, 0, 1));
    ASSERT_EQ(std::numeric_limits<GlobalIndexType>::max(), dof_map->getGlobalIndex(l_node0, 1, 0));

    MeshLib::Location l_node1(mesh->getID(), MeshLib::MeshItemType::Node, 1);
    ASSERT_EQ(1, dof_map->getGlobalIndex(l_node1, 0, 0));
    ASSERT_EQ(11, dof_map->getGlobalIndex(l_node1, 0, 1));
    ASSERT_EQ(20, dof_map->getGlobalIndex(l_node1, 1, 0));

    auto ele0_c0_indices = (*dof_map)(0, 0);
    ASSERT_EQ(2u, ele0_c0_indices.rows.size());
    auto ele0_c2_indices = (*dof_map)(0, 2);
    ASSERT_EQ(0u, ele0_c2_indices.rows.size());

    auto ele1_c2_indices = (*dof_map)(1, 2);
    ASSERT_EQ(2u, ele1_c2_indices.rows.size());
}


#ifndef USE_PETSC
TEST_F(NumLibLocalToGlobalIndexMapTest, MultipleVariablesMultipleComponentsHeterogeneousWithinElement)
#else
TEST_F(NumLibLocalToGlobalIndexMapTest, DISABLED_MultipleVariablesMultipleComponentsHeterogeneousWithinElement)
#endif
{
    // test 2 variables
    // - 1st variable with 2 components for all nodes, elements
    // - 2nd variable with 1 component for 1st node of element id 1
    std::vector<MeshLib::Node*> var2_nodes{const_cast<MeshLib::Node*>(mesh->getNode(1))};
    auto var2_subset =
        std::make_unique<MeshLib::MeshSubset>(*mesh, &var2_nodes);
    components.emplace_back(var2_subset.get());

    std::vector<int> vec_var_n_components{2, 1};
    std::vector<std::vector<MeshLib::Element*>const*> vec_var_elements;
    vec_var_elements.push_back(&mesh->getElements());
    std::vector<MeshLib::Element*> var2_elements{const_cast<MeshLib::Element*>(mesh->getElement(1))};
    vec_var_elements.push_back(&var2_elements);

    dof_map = std::make_unique<NumLib::LocalToGlobalIndexMap>(
        std::move(components),
        vec_var_n_components,
        vec_var_elements,
        NumLib::ComponentOrder::BY_COMPONENT);

    ASSERT_EQ(21u, dof_map->dofSizeWithGhosts());
    ASSERT_EQ(3, dof_map->getNumberOfComponents());
    ASSERT_EQ(2u, dof_map->getNumberOfVariables());
    ASSERT_EQ(2, dof_map->getNumberOfVariableComponents(0));
    ASSERT_EQ(1, dof_map->getNumberOfVariableComponents(1));

    MeshLib::Location l_node0(mesh->getID(), MeshLib::MeshItemType::Node, 0);

    ASSERT_EQ(0, dof_map->getGlobalIndex(l_node0, 0, 0));
    ASSERT_EQ(10, dof_map->getGlobalIndex(l_node0, 0, 1));
    ASSERT_EQ(std::numeric_limits<GlobalIndexType>::max(), dof_map->getGlobalIndex(l_node0, 1, 0));

    MeshLib::Location l_node1(mesh->getID(), MeshLib::MeshItemType::Node, 1);
    ASSERT_EQ(1, dof_map->getGlobalIndex(l_node1, 0, 0));
    ASSERT_EQ(11, dof_map->getGlobalIndex(l_node1, 0, 1));
    ASSERT_EQ(20, dof_map->getGlobalIndex(l_node1, 1, 0));

    auto ele0_c0_indices = (*dof_map)(0, 0);
    ASSERT_EQ(2u, ele0_c0_indices.rows.size());
    auto ele0_c2_indices = (*dof_map)(0, 2);
    ASSERT_EQ(0u, ele0_c2_indices.rows.size());

    auto ele1_c2_indices = (*dof_map)(1, 2);
    ASSERT_EQ(1u, ele1_c2_indices.rows.size());
    ASSERT_EQ(20u, ele1_c2_indices.rows[0]);
}
