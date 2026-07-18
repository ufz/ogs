// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <Eigen/Core>
#include <array>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Properties.h"
#include "MeshLib/PropertyVector.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"

class MeshLibProperties : public ::testing::Test
{
public:
    MeshLibProperties()
    {
        mesh =
            MeshToolsLib::MeshGenerator::generateRegularHexMesh(1.0, mesh_size);
    }

    ~MeshLibProperties() override { delete mesh; }

    /// Creates a single-component cell property of n_items values filled with
    /// first_value, first_value + 1, ... .
    MeshLib::PropertyVector<double>* createFilledCellProperty(
        std::string_view const name,
        std::size_t const n_items,
        double const first_value = 0.0)
    {
        auto* const pv = mesh->getProperties().createNewPropertyVector<double>(
            name, MeshLib::MeshItemType::Cell, n_items, 1);
        std::iota(pv->begin(), pv->end(), first_value);
        return pv;
    }

    static std::size_t const mesh_size = 5;
    MeshLib::Mesh* mesh{nullptr};
};
std::size_t const MeshLibProperties::mesh_size;

TEST_F(MeshLibProperties, PropertyVectorTestMetaData)
{
    ASSERT_TRUE(mesh != nullptr);

    std::string const prop_name("TestProperty");
    auto const* const p = mesh->getProperties().createNewPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::Cell);

    ASSERT_EQ(0u, p->getPropertyName().compare(prop_name));
    ASSERT_EQ(MeshLib::MeshItemType::Cell, p->getMeshItemType());
    ASSERT_EQ(1u, p->getNumberOfGlobalComponents());
    ASSERT_EQ(0u, p->size());
}

TEST_F(MeshLibProperties, PropertyVectorTestIntegrationPoint)
{
    ASSERT_TRUE(mesh != nullptr);

    int const n_integration_points = 4;
    int const vector_length = 6;
    std::string const prop_name("ip_field");

    auto* const p_ptr = mesh->getProperties().createNewPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::IntegrationPoint, vector_length);
    ASSERT_TRUE(p_ptr != nullptr);

    auto& p = *p_ptr;
    ASSERT_EQ(p.getPropertyName(), prop_name);
    ASSERT_EQ(MeshLib::MeshItemType::IntegrationPoint, p.getMeshItemType());
    ASSERT_EQ(vector_length, p.getNumberOfGlobalComponents());
    ASSERT_EQ(0u, p.size());

    // Fill the property vector with double data in following pattern:
    // for element's id xx and integration point number yy the value is xx.yy.
    p.resize(mesh->getNumberOfElements() * n_integration_points);

    std::vector<std::size_t> offsets(mesh->getNumberOfElements());
    std::size_t offset = 0;  // last position in the property vector
    for (auto const& e : mesh->getElements())
    {
        offsets[e->getID()] = offset;
        for (int ip = 0; ip < n_integration_points; ++ip)
        {
            p[offset + ip] = e->getID() + ip * 0.01;
        }
        offset += n_integration_points;
    }

    // Check the values at each offset.

    // Last element is checked after the for-loop.
    std::size_t i = 0;
    for (; i < offsets.size() - 1; ++i)
    {
        std::size_t const size = offsets[i + 1] - offsets[i];
        ASSERT_EQ(n_integration_points, size);
        for (int ip = 0; ip < n_integration_points; ++ip)
        {
            ASSERT_EQ(i + ip * 0.01, p[offsets[i] + ip]);
        }
    }
    {  // Last element
        std::size_t const size = p.size() - offsets[i];
        ASSERT_EQ(n_integration_points, size);
        for (int ip = 0; ip < n_integration_points; ++ip)
        {
            ASSERT_EQ(i + ip * 0.01, p[offsets[i] + ip]);
        }
    }
}

TEST_F(MeshLibProperties, AddDoubleProperties)
{
    ASSERT_TRUE(mesh != nullptr);
    const std::size_t size(mesh_size * mesh_size * mesh_size);

    std::string_view prop_name("TestProperty");
    auto* const double_properties =
        createFilledCellProperty(prop_name, size, 1.0);
    ASSERT_EQ(size, double_properties->size());

    for (std::size_t k(0); k < size; k++)
    {
        ASSERT_EQ(static_cast<double>(k + 1), (*double_properties)[k]);
    }

    ASSERT_TRUE(mesh->getProperties().existsPropertyVector<double>(prop_name));
    auto* const double_properties_cpy =
        mesh->getProperties().getPropertyVector<double>(prop_name);
    for (std::size_t k(0); k < size; k++)
    {
        ASSERT_EQ((*double_properties)[k], (*double_properties_cpy)[k]);
    }

    mesh->getProperties().removePropertyVector(prop_name);
    ASSERT_FALSE(mesh->getProperties().existsPropertyVector<double>(prop_name));
}

TEST_F(MeshLibProperties, AddVariousDifferentProperties)
{
    ASSERT_TRUE(mesh != nullptr);

    // *** add a 1st property ***
    std::string const& prop_name_2("ItemwiseMatrixProperties");
    // check if the property is already assigned to the mesh
    ASSERT_FALSE(mesh->getProperties().hasPropertyVector(prop_name_2));
    const std::size_t n_items_2(mesh_size * mesh_size * mesh_size);
    auto* const array_properties =
        mesh->getProperties().createNewPropertyVector<std::array<float, 9>>(
            prop_name_2, MeshLib::MeshItemType::Cell, n_items_2, 1);

    ASSERT_EQ(array_properties->size(), n_items_2);

    // initialize the property values
    for (std::size_t i(0); i < n_items_2; i++)
    {
        // init property value
        for (std::size_t k(0); k < (*array_properties)[i].size(); k++)
        {
            (*array_properties)[i][k] = static_cast<float>(i + k);
        }
    }

    EXPECT_EQ(9, (*array_properties)[0].size());

    // the mesh should have the property assigned to cells
    ASSERT_TRUE(mesh->getProperties().hasPropertyVector(prop_name_2));

    // fetch the vector in order to compare the content
    auto const* const array_properties_cpy =
        mesh->getProperties().getPropertyVector<std::array<float, 9>>(
            prop_name_2);
    ASSERT_FALSE(!array_properties_cpy);

    // compare the values/matrices
    for (std::size_t k(0); k < n_items_2; k++)
    {
        for (std::size_t j(0); j < (*array_properties)[k].size(); j++)
        {
            ASSERT_EQ((*array_properties)[k][j], (*array_properties_cpy)[k][j]);
        }
    }

    // *** add a 2nd property ***
    std::string const& prop_name_3("ItemwiseEigenMatrixProperties");
    // check if the property is already assigned to the mesh
    ASSERT_FALSE(mesh->getProperties().hasPropertyVector(prop_name_3));
    auto* const matrix_properties =
        mesh->getProperties()
            .createNewPropertyVector<
                Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(
                prop_name_3, MeshLib::MeshItemType::Cell);
    // init property values
    for (auto it = matrix_properties->begin(); it != matrix_properties->end();
         it++)
    {
        for (int r(0); r < it->rows(); r++)
        {
            for (int c(0); c < it->cols(); c++)
            {
                (*it)(r, c) = static_cast<double>(
                    std::distance(matrix_properties->begin(), it) +
                    r * it->cols() + c + 1);
            }
        }
    }

    // the mesh should have the property assigned to cells
    ASSERT_TRUE(mesh->getProperties().hasPropertyVector(prop_name_3));

    // fetch the vector in order to compare the content
    auto const* const matrix_properties_cpy =
        mesh->getProperties()
            .getPropertyVector<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(
                prop_name_3);
    ASSERT_FALSE(!matrix_properties_cpy);

    // compare the values/matrices
    auto it_cpy = matrix_properties_cpy->begin();
    for (auto it = matrix_properties->begin(); it != matrix_properties->end();
         it++, it_cpy++)
    {
        for (int r(0); r < it->rows(); r++)
        {
            for (int c(0); c < it->cols(); c++)
            {
                ASSERT_EQ((*it)(r, c), (*it_cpy)(r, c));
            }
        }
    }
}

TEST_F(MeshLibProperties, CopyConstructor)
{
    ASSERT_TRUE(mesh != nullptr);
    std::string const prop_name("CellProperty");
    const std::size_t n_items(mesh_size * mesh_size * mesh_size);

    // obtain PropertyVector data structure and fill it
    auto* const properties(createFilledCellProperty(prop_name, n_items, 1.0));

    // create a copy from the original Properties object
    MeshLib::Properties properties_copy(mesh->getProperties());
    // check if the Properties have a PropertyVector with the correct name
    ASSERT_TRUE(properties_copy.hasPropertyVector(prop_name));
    // fetch the PropertyVector from the copy of the Properties object
    auto const* const properties_cpy(
        properties_copy.getPropertyVector<double>(prop_name));
    ASSERT_FALSE(!properties_cpy);
    ASSERT_NE(properties, properties_cpy);  // distinct objects (deep copy)

    // check if the values in the PropertyVector of the copy of the Properties
    // are the same
    for (std::size_t k(0); k < n_items; k++)
    {
        EXPECT_EQ((*properties)[k], (*properties_cpy)[k]);
    }
}

TEST_F(MeshLibProperties, AddDoublePropertiesTupleSize2)
{
    ASSERT_TRUE(mesh != nullptr);
    const std::size_t number_of_tuples(mesh_size * mesh_size * mesh_size);

    std::string const prop_name("TestProperty");
    auto* const pv = mesh->getProperties().createNewPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::Cell, 2);
    // PropertyVector should be created in a correct way
    ASSERT_NE(nullptr, pv);

    ASSERT_EQ(0u, pv->getPropertyName().compare(prop_name));
    ASSERT_EQ(MeshLib::MeshItemType::Cell, pv->getMeshItemType());
    ASSERT_EQ(2u, pv->getNumberOfGlobalComponents());
    ASSERT_EQ(0u, pv->getNumberOfTuples());
    ASSERT_EQ(0u, pv->size());

    // push some values (2 tuples) into the vector
    for (std::size_t k(0); k < number_of_tuples; k++)
    {
        pv->push_back(static_cast<double>(k));
        pv->push_back(static_cast<double>(k));
    }
    // check the number of tuples
    ASSERT_EQ(number_of_tuples, pv->getNumberOfTuples());
    ASSERT_EQ(pv->getNumberOfTuples() * pv->getNumberOfGlobalComponents(),
              pv->size());
    // check the values
    for (std::size_t k(0); k < number_of_tuples; k++)
    {
        ASSERT_EQ(static_cast<double>(k), (*pv)[2 * k]);
        ASSERT_EQ(static_cast<double>(k), (*pv)[2 * k + 1]);
    }
}

// creation never overwrites; a second creation under the same name returns
// nullptr (and leaves the existing vector untouched).
TEST_F(MeshLibProperties, CreateWithDuplicateNameReturnsNullptr)
{
    ASSERT_TRUE(mesh != nullptr);
    std::string const prop_name("Duplicate");

    auto* const first = mesh->getProperties().createNewPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::Cell);
    ASSERT_NE(nullptr, first);

    // Same name, same type.
    auto* const second = mesh->getProperties().createNewPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::Cell);
    ASSERT_EQ(nullptr, second);

    // Same name, different type is still rejected: the map key is the name
    // alone.
    auto* const third = mesh->getProperties().createNewPropertyVector<int>(
        prop_name, MeshLib::MeshItemType::Cell);
    ASSERT_EQ(nullptr, third);
}

// getComponent gives 2-D (tuple, component) access into the flat,
// tuple-contiguous layout.
TEST_F(MeshLibProperties, GetComponentTwoDimensionalAccess)
{
    ASSERT_TRUE(mesh != nullptr);
    std::size_t const n_tuples(4);
    int const n_components(3);

    auto* const pv = mesh->getProperties().createNewPropertyVector<double>(
        "TwoD", MeshLib::MeshItemType::Cell, n_tuples, n_components);
    ASSERT_NE(nullptr, pv);
    ASSERT_EQ(n_tuples, pv->getNumberOfTuples());

    // Fill the flat data with 0, 1, 2, ... .
    std::iota(pv->begin(), pv->end(), 0.0);
    for (std::size_t t = 0; t < n_tuples; ++t)
    {
        for (int c = 0; c < n_components; ++c)
        {
            // Flat layout: data_[t * n_components + c].
            ASSERT_EQ((*pv)[t * n_components + c], pv->getComponent(t, c));
        }
    }

    // Mutating through getComponent is visible in the flat view.
    pv->getComponent(2, 1) = 42.0;
    ASSERT_EQ(42.0, (*pv)[2 * n_components + 1]);
}

// getPropertyVector is fail-fast - it OGS_FATALs (throws in tests) on a missing
// name, wrong type, wrong item type, or wrong component count.
TEST_F(MeshLibProperties, GetPropertyVectorMismatchThrows)
{
    ASSERT_TRUE(mesh != nullptr);
    std::string const prop_name("Typed");
    mesh->getProperties().createNewPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::Cell, 5, 2);

    auto& props = mesh->getProperties();

    // Missing name.
    EXPECT_THROW(props.getPropertyVector<double>("does_not_exist"),
                 std::runtime_error);
    // Wrong value type.
    EXPECT_THROW(props.getPropertyVector<int>(prop_name), std::runtime_error);
    // Wrong item type.
    EXPECT_THROW(props.getPropertyVector<double>(
                     prop_name, MeshLib::MeshItemType::Node, 2),
                 std::runtime_error);
    // Wrong component count.
    EXPECT_THROW(props.getPropertyVector<double>(
                     prop_name, MeshLib::MeshItemType::Cell, 1),
                 std::runtime_error);
    // Exact match succeeds.
    EXPECT_NO_THROW(props.getPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::Cell, 2));
}

// existsPropertyVector<T>(name, item_type, n_components) is a non-throwing
// probe that checks type, item type, and component count.
TEST_F(MeshLibProperties, ExistsPropertyVectorWithFullSignature)
{
    ASSERT_TRUE(mesh != nullptr);
    std::string const prop_name("Probe");
    mesh->getProperties().createNewPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::Cell, 5, 2);

    auto const& props = mesh->getProperties();

    EXPECT_TRUE(props.existsPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::Cell, 2));
    // Wrong type / item type / components each probe false without throwing.
    EXPECT_FALSE(props.existsPropertyVector<int>(
        prop_name, MeshLib::MeshItemType::Cell, 2));
    EXPECT_FALSE(props.existsPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::Node, 2));
    EXPECT_FALSE(props.existsPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::Cell, 1));
    EXPECT_FALSE(props.existsPropertyVector<double>(
        "missing", MeshLib::MeshItemType::Cell, 2));
}

// hasPropertyVector<T>(name, item_type) checks presence, type, and item type
// together.
TEST_F(MeshLibProperties, HasPropertyVectorWithItemType)
{
    ASSERT_TRUE(mesh != nullptr);
    std::string const prop_name("Field");
    mesh->getProperties().createNewPropertyVector<double>(
        prop_name, MeshLib::MeshItemType::Cell);

    auto const& props = mesh->getProperties();
    EXPECT_TRUE(props.hasPropertyVector<double>(prop_name,
                                                MeshLib::MeshItemType::Cell));
    EXPECT_FALSE(props.hasPropertyVector<double>(prop_name,
                                                 MeshLib::MeshItemType::Node));
    EXPECT_FALSE(
        props.hasPropertyVector<int>(prop_name, MeshLib::MeshItemType::Cell));
}

// name enumeration, globally and filtered by item type.
TEST_F(MeshLibProperties, EnumerateNamesAndSizes)
{
    ASSERT_TRUE(mesh != nullptr);
    auto& props = mesh->getProperties();
    props.createNewPropertyVector<double>("cell_a",
                                          MeshLib::MeshItemType::Cell);
    props.createNewPropertyVector<double>("cell_b",
                                          MeshLib::MeshItemType::Cell);
    props.createNewPropertyVector<int>("node_a", MeshLib::MeshItemType::Node);

    EXPECT_EQ(3u, props.size());
    EXPECT_EQ(2u, props.size(MeshLib::MeshItemType::Cell));
    EXPECT_EQ(1u, props.size(MeshLib::MeshItemType::Node));

    auto const all_names = props.getPropertyVectorNames();
    EXPECT_EQ(3u, all_names.size());

    auto const cell_names =
        props.getPropertyVectorNames(MeshLib::MeshItemType::Cell);
    ASSERT_EQ(2u, cell_names.size());
    EXPECT_EQ("cell_a", cell_names[0]);
    EXPECT_EQ("cell_b", cell_names[1]);
}

// assign(R&&) bulk-assigns from an arbitrary input range.
TEST_F(MeshLibProperties, AssignFromRange)
{
    ASSERT_TRUE(mesh != nullptr);
    auto* const pv = mesh->getProperties().createNewPropertyVector<double>(
        "Assigned", MeshLib::MeshItemType::Cell);
    ASSERT_TRUE(pv->empty());

    std::vector<double> const source{1.5, 2.5, 3.5, 4.5};
    pv->assign(source);

    ASSERT_EQ(source.size(), pv->size());
    for (std::size_t i = 0; i < source.size(); ++i)
    {
        EXPECT_EQ(source[i], (*pv)[i]);
    }

    ASSERT_EQ(std::ssize(source), pv->ssize());
}

// clear() empties, empty() reflects state.
TEST_F(MeshLibProperties, ClearAndEmpty)
{
    ASSERT_TRUE(mesh != nullptr);
    auto* const pv = mesh->getProperties().createNewPropertyVector<double>(
        "Clearable", MeshLib::MeshItemType::Cell, 3, 1);
    ASSERT_FALSE(pv->empty());
    ASSERT_EQ(3u, pv->size());

    pv->clear();
    EXPECT_TRUE(pv->empty());
    EXPECT_EQ(0u, pv->size());
}

// Properties copy assignment deep-copies the primary template, so mutating the
// original does not affect the copy.
TEST_F(MeshLibProperties, CopyAssignmentDeepCopiesPrimaryTemplate)
{
    ASSERT_TRUE(mesh != nullptr);
    std::string const prop_name("DeepCopy");
    auto* const original = createFilledCellProperty(prop_name, 4, 1.0);

    MeshLib::Properties copy;
    copy = mesh->getProperties();

    auto* const copied = copy.getPropertyVector<double>(prop_name);
    ASSERT_NE(original, copied);  // distinct objects
    ASSERT_EQ(original->size(), copied->size());

    // Mutate the original; the copy must be unaffected.
    (*original)[0] = 999.0;
    EXPECT_EQ(1.0, (*copied)[0]);
    EXPECT_EQ(999.0, (*original)[0]);
}

// excludeCopyProperties(elem_ids, node_ids) clones Cell/Node fields dropping
// the requested tuples (single-component case).
TEST_F(MeshLibProperties, ExcludeCopyPropertiesByIds)
{
    ASSERT_TRUE(mesh != nullptr);
    createFilledCellProperty("cell", 5);  // 0,1,2,3,4

    // Ascending order is required by excludeObjectCopy.
    std::vector<std::size_t> const exclude_elem_ids{1, 3};
    std::vector<std::size_t> const exclude_node_ids{};

    MeshLib::Properties const reduced =
        mesh->getProperties().excludeCopyProperties(exclude_elem_ids,
                                                    exclude_node_ids);

    auto const* const reduced_cell = reduced.getPropertyVector<double>("cell");
    ASSERT_EQ(3u, reduced_cell->size());  // 5 - 2 excluded
    EXPECT_EQ(0.0, (*reduced_cell)[0]);
    EXPECT_EQ(2.0, (*reduced_cell)[1]);
    EXPECT_EQ(4.0, (*reduced_cell)[2]);
}

// excludeCopyProperties(item_types) copies every field except those whose item
// type is excluded (no tuple dropping).
TEST_F(MeshLibProperties, ExcludeCopyPropertiesByItemType)
{
    ASSERT_TRUE(mesh != nullptr);
    auto& props = mesh->getProperties();
    props.createNewPropertyVector<double>("cell", MeshLib::MeshItemType::Cell,
                                          3, 1);
    props.createNewPropertyVector<int>("node", MeshLib::MeshItemType::Node, 3,
                                       1);

    MeshLib::Properties const kept = props.excludeCopyProperties(
        std::vector<MeshLib::MeshItemType>{MeshLib::MeshItemType::Node});

    EXPECT_TRUE(kept.hasPropertyVector("cell"));
    EXPECT_FALSE(kept.hasPropertyVector("node"));
    EXPECT_EQ(1u, kept.size());
}

// applyToPropertyVectors visits each property with its resolved value type
// exactly once.
TEST_F(MeshLibProperties, ApplyToPropertyVectorsVisitsEachOnce)
{
    ASSERT_TRUE(mesh != nullptr);
    auto& props = mesh->getProperties();
    props.createNewPropertyVector<double>("d", MeshLib::MeshItemType::Cell, 2,
                                          1);
    props.createNewPropertyVector<int>("i", MeshLib::MeshItemType::Cell, 2, 1);

    int visited = 0;
    applyToPropertyVectors(
        props,
        [&visited](auto type, MeshLib::PropertyVectorBase* const property)
        {
            auto* const p =
                dynamic_cast<MeshLib::PropertyVector<decltype(type)>*>(
                    property);
            if (p == nullptr)
            {
                return false;
            }
            ++visited;
            return true;
        });

    EXPECT_EQ(2, visited);
}
