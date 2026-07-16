// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "MathLib/Point3d.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Utils/IntegrationPointMetaData.h"
#include "MeshLib/Utils/IntegrationPointWriter.h"
#include "MeshLib/Utils/getOrCreateMeshProperty.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"
#include "MeshToolsLib/OverwriteMeshFieldDataByMaterialIDs.h"
#include "ParameterLib/ConstantParameter.h"
#include "ParameterLib/SpatialPosition.h"

namespace
{
// Helper function to create a simple line mesh with material IDs
std::unique_ptr<MeshLib::Mesh> createLineMeshWithMaterialIDs(
    std::size_t num_elements, std::vector<int> const& material_ids)
{
    auto mesh = std::unique_ptr<MeshLib::Mesh>(
        MeshToolsLib::MeshGenerator::generateLineMesh(
            static_cast<double>(num_elements), num_elements, MathLib::ORIGIN,
            "test_mesh"));

    auto* mat_ids = MeshLib::getOrCreateMeshProperty<int>(
        *mesh, "MaterialIDs", MeshLib::MeshItemType::Cell, 1);

    for (std::size_t i = 0; i < material_ids.size() && i < mat_ids->size(); ++i)
    {
        (*mat_ids)[i] = material_ids[i];
    }

    return mesh;
}

// Helper to get element IDs for a given material ID
std::vector<std::size_t> getElementIDsForMaterial(MeshLib::Mesh const& mesh,
                                                  int material_id)
{
    auto const* mat_ids = mesh.getProperties().getPropertyVector<int>(
        "MaterialIDs", MeshLib::MeshItemType::Cell, 1);

    std::vector<std::size_t> element_ids;
    if (mat_ids)
    {
        for (std::size_t i = 0; i < mat_ids->size(); ++i)
        {
            if ((*mat_ids)[i] == material_id)
            {
                element_ids.push_back(i);
            }
        }
    }
    return element_ids;
}
}  // namespace

// Test fixture for cell data tests
class OverwriteCellDataTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create a line mesh with 10 elements
        mesh =
            createLineMeshWithMaterialIDs(10, {0, 0, 0, 0, 0, 1, 1, 1, 1, 1});

        // Create a cell property to test overwriting
        property_name = "TestProperty";
        auto* prop = MeshLib::getOrCreateMeshProperty<double>(
            *mesh, property_name, MeshLib::MeshItemType::Cell, 1);

        // Initialize with distinct values
        for (std::size_t i = 0; i < prop->size(); ++i)
        {
            (*prop)[i] = 100.0 + static_cast<double>(i);
        }
    }

    std::unique_ptr<MeshLib::Mesh> mesh;
    std::string property_name;
};

// Test fixture for node data tests
class OverwriteNodeDataTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create a line mesh with 10 elements (11 nodes)
        mesh =
            createLineMeshWithMaterialIDs(10, {0, 0, 0, 0, 0, 1, 1, 1, 1, 1});

        // Create a node property to test overwriting
        property_name = "TestProperty";
        auto* prop = MeshLib::getOrCreateMeshProperty<double>(
            *mesh, property_name, MeshLib::MeshItemType::Node, 1);

        // Initialize with distinct values
        for (std::size_t i = 0; i < prop->size(); ++i)
        {
            (*prop)[i] = 1000.0 + static_cast<double>(i);
        }
    }

    std::unique_ptr<MeshLib::Mesh> mesh;
    std::string property_name;
};

// Test overwriting cell data with constant value (no parameter)
TEST_F(OverwriteCellDataTest, OverwriteWithConstantValue)
{
    auto const element_ids = getElementIDsForMaterial(*mesh, 0);
    ASSERT_EQ(5u, element_ids.size());

    // Setup initial condition data
    MeshToolsLib::InitialConditionDataSet initial_condition;
    initial_condition.variable_name = property_name;
    initial_condition.element_ids_for_selected_materials = element_ids;
    initial_condition.mesh = mesh.get();
    initial_condition.parameter = nullptr;
    initial_condition.mesh_item_type = MeshLib::MeshItemType::Cell;

    // Perform the overwrite using the public API
    std::vector<MeshToolsLib::InitialConditionDataSet> data;
    data.push_back(std::move(initial_condition));
    MeshToolsLib::overwriteMeshFieldDataByMaterialIDs(data);

    // Verify that elements 0-4 now have their original values (since each
    // element copies its own value from the property copy)
    auto* updated_prop = MeshLib::getOrCreateMeshProperty<double>(
        *mesh, property_name, MeshLib::MeshItemType::Cell, 1);

    for (std::size_t i = 0; i < 5; ++i)
    {
        ASSERT_EQ(100.0 + i, (*updated_prop)[i])
            << "Element " << i << " should have its original value";
    }

    // Elements 5-9 should remain unchanged
    for (std::size_t i = 5; i < 10; ++i)
    {
        ASSERT_EQ(100.0 + i, (*updated_prop)[i])
            << "Element " << i << " should remain unchanged";
    }
}

// Test overwriting cell data with parameter value
TEST_F(OverwriteCellDataTest, OverwriteWithParameterValue)
{
    auto const element_ids = getElementIDsForMaterial(*mesh, 0);
    ASSERT_EQ(5u, element_ids.size());

    // Create a constant parameter that returns 42.0
    ParameterLib::ConstantParameter<double> param{"test_param", 42.0};

    // Setup initial condition data with parameter
    MeshToolsLib::InitialConditionDataSet initial_condition;
    initial_condition.variable_name = property_name;
    initial_condition.element_ids_for_selected_materials = element_ids;
    initial_condition.mesh = mesh.get();
    initial_condition.parameter = &param;
    initial_condition.mesh_item_type = MeshLib::MeshItemType::Cell;

    // Perform the overwrite using the public API
    std::vector<MeshToolsLib::InitialConditionDataSet> data;
    data.push_back(std::move(initial_condition));
    MeshToolsLib::overwriteMeshFieldDataByMaterialIDs(data);

    // Verify that elements 0-4 now have the parameter value (42.0)
    auto* prop = MeshLib::getOrCreateMeshProperty<double>(
        *mesh, property_name, MeshLib::MeshItemType::Cell, 1);

    for (std::size_t i = 0; i < 5; ++i)
    {
        ASSERT_EQ(42.0, (*prop)[i])
            << "Element " << i << " should have value 42.0 from parameter";
    }

    // Elements 5-9 should remain unchanged
    for (std::size_t i = 5; i < 10; ++i)
    {
        ASSERT_EQ(100.0 + i, (*prop)[i])
            << "Element " << i << " should remain unchanged";
    }
}

// Test overwriting node data with constant value (no parameter)
TEST_F(OverwriteNodeDataTest, OverwriteNodeDataWithConstantValue)
{
    auto const element_ids = getElementIDsForMaterial(*mesh, 0);
    ASSERT_EQ(5u, element_ids.size());

    // Setup initial condition data
    MeshToolsLib::InitialConditionDataSet initial_condition;
    initial_condition.variable_name = property_name;
    initial_condition.element_ids_for_selected_materials = element_ids;
    initial_condition.mesh = mesh.get();
    initial_condition.parameter = nullptr;
    initial_condition.mesh_item_type = MeshLib::MeshItemType::Node;

    // Perform the overwrite using the public API
    std::vector<MeshToolsLib::InitialConditionDataSet> data;
    data.push_back(std::move(initial_condition));
    MeshToolsLib::overwriteMeshFieldDataByMaterialIDs(data);

    // Verify that nodes belonging to elements 0-4 have been updated
    // Since each node copies its own value from the property copy,
    // nodes 0-5 should have their original values (1000.0, 1001.0, ..., 1005.0)
    auto* updated_prop = MeshLib::getOrCreateMeshProperty<double>(
        *mesh, property_name, MeshLib::MeshItemType::Node, 1);

    // Nodes 0-5 should have their original values
    for (std::size_t i = 0; i <= 5; ++i)
    {
        ASSERT_EQ(1000.0 + i, (*updated_prop)[i])
            << "Node " << i << " should have its original value";
    }
}

// Test overwriting node data with parameter value
TEST_F(OverwriteNodeDataTest, OverwriteNodeDataWithParameterValue)
{
    auto const element_ids = getElementIDsForMaterial(*mesh, 0);
    ASSERT_EQ(5u, element_ids.size());

    // Create a constant parameter that returns 77.0
    ParameterLib::ConstantParameter<double> param{"test_param", 77.0};

    // Setup initial condition data with parameter
    MeshToolsLib::InitialConditionDataSet initial_condition;
    initial_condition.variable_name = property_name;
    initial_condition.element_ids_for_selected_materials = element_ids;
    initial_condition.mesh = mesh.get();
    initial_condition.parameter = &param;
    initial_condition.mesh_item_type = MeshLib::MeshItemType::Node;

    // Perform the overwrite using the public API
    std::vector<MeshToolsLib::InitialConditionDataSet> data;
    data.push_back(std::move(initial_condition));
    MeshToolsLib::overwriteMeshFieldDataByMaterialIDs(data);

    // Verify that nodes belonging to elements 0-4 have been updated
    auto* prop = MeshLib::getOrCreateMeshProperty<double>(
        *mesh, property_name, MeshLib::MeshItemType::Node, 1);

    // Nodes 0-5 should have the parameter value (element 0-4 cover nodes 0-5)
    for (std::size_t i = 0; i <= 5; ++i)
    {
        ASSERT_EQ(77.0, (*prop)[i])
            << "Node " << i << " should have value 77.0 from parameter";
    }
}

// Test with multi-component property
TEST_F(OverwriteCellDataTest, OverwriteMultiComponentProperty)
{
    auto const element_ids = getElementIDsForMaterial(*mesh, 0);
    ASSERT_EQ(5u, element_ids.size());

    // Create a multi-component property (3 components)
    std::string multi_prop_name = "MultiComponentProperty";
    auto* prop = MeshLib::getOrCreateMeshProperty<double>(
        *mesh, multi_prop_name, MeshLib::MeshItemType::Cell, 3);

    // Initialize with distinct values for each component
    for (std::size_t i = 0; i < prop->size(); ++i)
    {
        (*prop)[i] = 1000.0 + static_cast<double>(i);
    }

    // Setup initial condition data
    MeshToolsLib::InitialConditionDataSet initial_condition;
    initial_condition.variable_name = multi_prop_name;
    initial_condition.element_ids_for_selected_materials = element_ids;
    initial_condition.mesh = mesh.get();
    initial_condition.parameter = nullptr;
    initial_condition.mesh_item_type = MeshLib::MeshItemType::Cell;

    // Perform the overwrite using the public API
    std::vector<MeshToolsLib::InitialConditionDataSet> data;
    data.push_back(std::move(initial_condition));
    MeshToolsLib::overwriteMeshFieldDataByMaterialIDs(data);

    // Verify that elements 0-4 have their original values
    // (each element copies its own values from the property copy)
    auto* updated_prop = MeshLib::getOrCreateMeshProperty<double>(
        *mesh, multi_prop_name, MeshLib::MeshItemType::Cell, 3);

    for (std::size_t i = 0; i < 5; ++i)
    {
        ASSERT_EQ(1000.0 + i * 3 + 0, (*updated_prop)[i * 3 + 0])
            << "Element " << i << " component 0";
        ASSERT_EQ(1000.0 + i * 3 + 1, (*updated_prop)[i * 3 + 1])
            << "Element " << i << " component 1";
        ASSERT_EQ(1000.0 + i * 3 + 2, (*updated_prop)[i * 3 + 2])
            << "Element " << i << " component 2";
    }
}

// Test that non-selected elements remain unchanged
TEST_F(OverwriteCellDataTest, NonSelectedElementsUnchanged)
{
    auto const element_ids_0 = getElementIDsForMaterial(*mesh, 0);
    auto const element_ids_1 = getElementIDsForMaterial(*mesh, 1);

    ASSERT_EQ(5u, element_ids_0.size());
    ASSERT_EQ(5u, element_ids_1.size());

    // Setup initial condition data for material 0 only
    MeshToolsLib::InitialConditionDataSet initial_condition;
    initial_condition.variable_name = property_name;
    initial_condition.element_ids_for_selected_materials = element_ids_0;
    initial_condition.mesh = mesh.get();
    initial_condition.parameter = nullptr;
    initial_condition.mesh_item_type = MeshLib::MeshItemType::Cell;

    // Perform the overwrite using the public API
    std::vector<MeshToolsLib::InitialConditionDataSet> data;
    data.push_back(std::move(initial_condition));
    MeshToolsLib::overwriteMeshFieldDataByMaterialIDs(data);

    // Verify elements 0-4 have their original values (each element copies its
    // own value)
    auto* updated_prop = MeshLib::getOrCreateMeshProperty<double>(
        *mesh, property_name, MeshLib::MeshItemType::Cell, 1);
    for (std::size_t i = 0; i < 5; ++i)
    {
        ASSERT_EQ(100.0 + i, (*updated_prop)[i])
            << "Element " << i << " should have its original value";
    }

    // Verify elements 5-9 (material 1) remain unchanged
    for (std::size_t i = 5; i < 10; ++i)
    {
        ASSERT_EQ(100.0 + i, (*updated_prop)[i])
            << "Element " << i << " should remain unchanged";
    }
}

// Test with empty element ID list
TEST_F(OverwriteCellDataTest, EmptyElementIDList)
{
    // Setup initial condition data with empty element list
    MeshToolsLib::InitialConditionDataSet initial_condition;
    initial_condition.variable_name = property_name;
    initial_condition.element_ids_for_selected_materials = {};
    initial_condition.mesh = mesh.get();
    initial_condition.parameter = nullptr;
    initial_condition.mesh_item_type = MeshLib::MeshItemType::Cell;

    // Perform the overwrite using the public API (should do nothing)
    std::vector<MeshToolsLib::InitialConditionDataSet> data;
    data.push_back(std::move(initial_condition));
    MeshToolsLib::overwriteMeshFieldDataByMaterialIDs(data);

    // Verify all elements remain unchanged
    auto* updated_prop = MeshLib::getOrCreateMeshProperty<double>(
        *mesh, property_name, MeshLib::MeshItemType::Cell, 1);
    for (std::size_t i = 0; i < 10; ++i)
    {
        ASSERT_EQ(100.0 + i, (*updated_prop)[i]);
    }
}

// Helper function to set up integration point metadata for a property
void setupIntegrationPointMetaData(MeshLib::Mesh& mesh,
                                   std::string const& property_name,
                                   int n_components,
                                   int integration_order = 2)
{
    // Create the IntegrationPointMetaData as a JSON string
    MeshLib::IntegrationPointMetaDataSingleField single_field{
        property_name, n_components, integration_order};
    std::vector<MeshLib::IntegrationPointMetaDataSingleField> fields{
        single_field};
    MeshLib::IntegrationPointMetaData ip_meta_data{fields};

    // Store it as a char property vector
    auto& dictionary = *MeshLib::getOrCreateMeshProperty<char>(
        mesh, "IntegrationPointMetaData",
        MeshLib::MeshItemType::IntegrationPoint, 1);
    dictionary.clear();
    std::string const json_string = ip_meta_data.toJsonString();
    std::copy(json_string.begin(), json_string.end(),
              std::back_inserter(dictionary));
}

// Test fixture for integration point data tests
class OverwriteIntegrationPointDataTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Create a hex mesh (2x2x2 = 8 elements)
        mesh = std::unique_ptr<MeshLib::Mesh>(
            MeshToolsLib::MeshGenerator::generateRegularHexMesh(1.0, 2));

        // Set up material IDs: first 4 elements are material 0, last 4 are
        // material 1
        auto* mat_ids = MeshLib::getOrCreateMeshProperty<int>(
            *mesh, "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
        for (std::size_t i = 0; i < 4; ++i)
        {
            (*mat_ids)[i] = 0;
        }
        for (std::size_t i = 4; i < mesh->getNumberOfElements(); ++i)
        {
            (*mat_ids)[i] = 1;
        }

        // Create an integration point property
        prop_name = "test_ip_property";
        ip_prop = mesh->getProperties().createNewPropertyVector<double>(
            prop_name, MeshLib::MeshItemType::IntegrationPoint, 1);

        // Resize to hold all integration points (8 elements * 8 IPs = 64
        // values) For hex elements with integration order 2, there are 2x2x2 =
        // 8 Gauss points
        std::size_t const n_elements = mesh->getNumberOfElements();
        std::size_t const n_ips_per_element = 8;
        ip_prop->resize(n_elements * n_ips_per_element);

        // Set up integration point metadata
        setupIntegrationPointMetaData(*mesh, prop_name, 1, 2);

        // Initialize with distinct values
        for (std::size_t i = 0; i < ip_prop->size(); ++i)
        {
            (*ip_prop)[i] = 1000.0 + static_cast<double>(i);
        }

        // Get element IDs for material 0
        element_ids.clear();
        for (std::size_t i = 0; i < 4; ++i)
        {
            element_ids.push_back(i);
        }
    }

    std::unique_ptr<MeshLib::Mesh> mesh;
    std::string prop_name;
    MeshLib::PropertyVector<double>* ip_prop = nullptr;
    std::vector<std::size_t> element_ids;
};

// Test for integration point data with constant value (copy mode)
TEST_F(OverwriteIntegrationPointDataTest, OverwriteWithConstantValue)
{
    // Setup initial condition data
    MeshToolsLib::InitialConditionDataSet initial_condition;
    initial_condition.variable_name = prop_name;
    initial_condition.element_ids_for_selected_materials = element_ids;
    initial_condition.mesh = mesh.get();
    initial_condition.parameter = nullptr;
    initial_condition.mesh_item_type = MeshLib::MeshItemType::IntegrationPoint;

    // Perform the overwrite using the public API
    std::vector<MeshToolsLib::InitialConditionDataSet> data;
    data.push_back(std::move(initial_condition));
    MeshToolsLib::overwriteMeshFieldDataByMaterialIDs(data);

    // Verify that integration points for elements 0-3 have their original
    // values (each IP copies its own value from the property copy)
    auto* updated_ip_prop = mesh->getProperties().getPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::IntegrationPoint, 1);

    // Each element has 8 IPs, so elements 0-3 have IPs 0-31
    for (std::size_t i = 0; i < 32; ++i)
    {
        ASSERT_EQ(1000.0 + i, (*updated_ip_prop)[i])
            << "IP " << i << " should have its original value";
    }

    // IPs for elements 4-7 should remain unchanged
    for (std::size_t i = 32; i < 64; ++i)
    {
        ASSERT_EQ(1000.0 + i, (*updated_ip_prop)[i])
            << "IP " << i << " should remain unchanged";
    }
}

// Test for integration point data with parameter
TEST_F(OverwriteIntegrationPointDataTest, OverwriteWithParameter)
{
    // Create a constant parameter
    ParameterLib::ConstantParameter<double> param{"test_param", 777.0};

    // Setup initial condition data with parameter
    MeshToolsLib::InitialConditionDataSet initial_condition;
    initial_condition.variable_name = prop_name;
    initial_condition.element_ids_for_selected_materials = element_ids;
    initial_condition.mesh = mesh.get();
    initial_condition.parameter = &param;
    initial_condition.mesh_item_type = MeshLib::MeshItemType::IntegrationPoint;

    // Perform the overwrite using the public API
    std::vector<MeshToolsLib::InitialConditionDataSet> data;
    data.push_back(std::move(initial_condition));
    MeshToolsLib::overwriteMeshFieldDataByMaterialIDs(data);

    // Verify that integration points for elements 0-3 have the parameter value
    auto* updated_ip_prop = mesh->getProperties().getPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::IntegrationPoint, 1);

    for (std::size_t i = 0; i < 32; ++i)
    {
        ASSERT_EQ(777.0, (*updated_ip_prop)[i])
            << "IP " << i << " should have value 777.0 from parameter";
    }

    // IPs for elements 4-7 should remain unchanged
    for (std::size_t i = 32; i < 64; ++i)
    {
        ASSERT_EQ(1000.0 + i, (*updated_ip_prop)[i])
            << "IP " << i << " should remain unchanged";
    }
}

// Test the main orchestration function
TEST(MeshToolsLib, OverwriteMeshFieldDataByMaterialIDsOrchestration)
{
    // Create a line mesh with 10 elements
    auto mesh = std::unique_ptr<MeshLib::Mesh>(
        MeshToolsLib::MeshGenerator::generateLineMesh(10.0, 10));

    // Set up material IDs
    auto* mat_ids = MeshLib::getOrCreateMeshProperty<int>(
        *mesh, "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
    for (std::size_t i = 0; i < 5; ++i)
    {
        (*mat_ids)[i] = 0;
    }
    for (std::size_t i = 5; i < 10; ++i)
    {
        (*mat_ids)[i] = 1;
    }

    // Get element IDs for material 0
    std::vector<std::size_t> element_ids_0 = {0, 1, 2, 3, 4};

    // Create a cell property
    std::string const prop_name = "TestProperty";
    auto* prop = MeshLib::getOrCreateMeshProperty<double>(
        *mesh, prop_name, MeshLib::MeshItemType::Cell, 1);
    for (std::size_t i = 0; i < prop->size(); ++i)
    {
        (*prop)[i] = 100.0 + static_cast<double>(i);
    }

    // Setup initial condition data for the main function
    MeshToolsLib::InitialConditionDataSet initial_condition;
    initial_condition.variable_name = prop_name;
    initial_condition.element_ids_for_selected_materials = element_ids_0;
    initial_condition.mesh = mesh.get();
    initial_condition.parameter = nullptr;
    initial_condition.mesh_item_type = MeshLib::MeshItemType::Cell;

    std::vector<MeshToolsLib::InitialConditionDataSet> data;
    data.push_back(std::move(initial_condition));

    // Call the orchestration function
    MeshToolsLib::overwriteMeshFieldDataByMaterialIDs(data);

    // The function copies the original property and then overwrites.
    // Each element copies its own value from the property_copy.
    // So elements 0-4 will keep their original values (100, 101, 102, 103, 104)
    auto* updated_prop = MeshLib::getOrCreateMeshProperty<double>(
        *mesh, prop_name, MeshLib::MeshItemType::Cell, 1);

    // Elements 0-4 should have their original values (100, 101, 102, 103, 104)
    // because each element copies its own value from the property copy
    for (std::size_t i = 0; i < 5; ++i)
    {
        ASSERT_EQ(100.0 + i, (*updated_prop)[i])
            << "Element " << i << " should have its original value";
    }

    // Elements 5-9 should remain unchanged
    for (std::size_t i = 5; i < 10; ++i)
    {
        ASSERT_EQ(100.0 + i, (*updated_prop)[i])
            << "Element " << i << " should remain unchanged";
    }
}
