/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <gtest/gtest.h>

#include <memory>
#include <range/v3/algorithm/fill.hpp>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/MeshEditing/AddLayerToMesh.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"
#include "MeshToolsLib/MeshInformation.h"
#include "MeshToolsLib/MeshQuality/MeshValidation.h"
#include "MeshToolsLib/MeshSurfaceExtraction.h"

namespace AddLayerValidation
{
// validates mesh. for line meshes, node order tests will fail because vertical
// elements are created during the extrusion. Therefore node order tests are
// switched off for line meshes.
void validate(MeshLib::Mesh const& mesh, bool testNodeOrder)
{
    int const reduce_tests = (testNodeOrder) ? 0 : 1;

    auto const nErrorFlags(
        static_cast<std::size_t>(ElementErrorFlag::MaxValue));
    ElementErrorFlag const flags[nErrorFlags] = {
        ElementErrorFlag::ZeroVolume, ElementErrorFlag::NonCoplanar,
        ElementErrorFlag::NonConvex, ElementErrorFlag::NodeOrder};
    std::vector<ElementErrorCode> const codes(
        MeshToolsLib::MeshValidation::testElementGeometry(mesh));
    for (auto code : codes)
    {
        for (std::size_t j = 0; j < nErrorFlags - reduce_tests; ++j)
        {
            ASSERT_FALSE(code[flags[j]]);
        }
    }
}

void testZCoords2D(MeshLib::Mesh const& input, MeshLib::Mesh const& output,
                   double height)
{
    std::size_t const nNodes(input.getNumberOfNodes());
    for (std::size_t i = 0; i < nNodes; ++i)
    {
        ASSERT_EQ((*input.getNode(i))[2], (*output.getNode(i))[2]);
        ASSERT_EQ((*input.getNode(i))[2] + height,
                  (*output.getNode(nNodes + i))[2]);
    }
}

void testZCoords3D(MeshLib::Mesh const& input, MeshLib::Mesh const& output,
                   double height)
{
    std::size_t const nNodes(input.getNumberOfNodes());
    for (std::size_t i = 0; i < nNodes; ++i)
    {
        ASSERT_EQ((*input.getNode(i))[2] + height, (*output.getNode(i))[2]);
    }
}
}  // namespace AddLayerValidation

TEST(MeshLib, AddTopLayerToLineMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshToolsLib::MeshGenerator::generateLineMesh(1.0, 5));
    double const height(1);
    constexpr bool const copy_material_ids = false;
    std::unique_ptr<MeshLib::Mesh> const result(MeshToolsLib::addLayerToMesh(
        *mesh, height, "mesh", true, copy_material_ids, std::nullopt));

    ASSERT_EQ(2 * mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(2 * mesh->getNumberOfElements(), result->getNumberOfElements());

    auto const& n_elems =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*result);
    ASSERT_EQ(5, n_elems.at(MeshLib::MeshElemType::LINE));
    ASSERT_EQ(5, n_elems.at(MeshLib::MeshElemType::QUAD));

    AddLayerValidation::testZCoords2D(*mesh, *result, height);
    AddLayerValidation::validate(*result, false);
}

TEST(MeshLib, AddBottomLayerToLineMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshToolsLib::MeshGenerator::generateLineMesh(1.0, 5));
    double const height(1);
    constexpr bool const copy_material_ids = false;
    std::unique_ptr<MeshLib::Mesh> const result(MeshToolsLib::addLayerToMesh(
        *mesh, height, "mesh", false, copy_material_ids, std::nullopt));

    ASSERT_EQ(2 * mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(2 * mesh->getNumberOfElements(), result->getNumberOfElements());

    auto const& n_elems =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*result);
    ASSERT_EQ(5, n_elems.at(MeshLib::MeshElemType::LINE));
    ASSERT_EQ(5, n_elems.at(MeshLib::MeshElemType::QUAD));

    AddLayerValidation::testZCoords2D(*mesh, *result, -1 * height);
    AddLayerValidation::validate(*result, false);
}

TEST(MeshLib, AddTopLayerToTriMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshToolsLib::MeshGenerator::generateRegularTriMesh(5, 5));
    std::string const& mat_name("MaterialIDs");
    auto* mats = mesh->getProperties().createNewPropertyVector<int>(
        mat_name, MeshLib::MeshItemType::Cell, mesh->getNumberOfElements(), 1);
    ASSERT_NE(nullptr, mats);
    if (mats)
    {
        ranges::fill(mats->begin(), mats->end(), 0);
    }
    auto* test = mesh->getProperties().createNewPropertyVector<double>(
        "test", MeshLib::MeshItemType::Cell, mesh->getNumberOfElements(), 1);
    ASSERT_NE(nullptr, test);
    if (test)
    {
        ranges::fill(*test, 0.1);
    }
    ASSERT_EQ(2, mesh->getProperties().size());
    double const height(1);
    constexpr bool copy_material_ids = false;
    std::unique_ptr<MeshLib::Mesh> const result(MeshToolsLib::addLayerToMesh(
        *mesh, height, "mesh", true, copy_material_ids, std::nullopt));

    ASSERT_EQ(2 * mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(2 * mesh->getNumberOfElements(), result->getNumberOfElements());

    auto const& n_elems =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*result);
    ASSERT_EQ(mesh->getNumberOfElements(),
              n_elems.at(MeshLib::MeshElemType::TRIANGLE));
    ASSERT_EQ(mesh->getNumberOfElements(),
              n_elems.at(MeshLib::MeshElemType::PRISM));

    ASSERT_EQ(1, result->getProperties().size());
    auto const* const new_mats =
        result->getProperties().getPropertyVector<int>(mat_name);
    ASSERT_EQ(result->getNumberOfElements(), new_mats->size());
    ASSERT_EQ(mesh->getNumberOfElements(),
              std::count(new_mats->cbegin(), new_mats->cend(), 0));
    ASSERT_EQ(mesh->getNumberOfElements(),
              std::count(new_mats->cbegin(), new_mats->cend(), 1));
    AddLayerValidation::testZCoords2D(*mesh, *result, height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddBottomLayerToTriMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshToolsLib::MeshGenerator::generateRegularTriMesh(5, 5));
    double const height(1);
    constexpr bool copy_material_ids = false;
    std::unique_ptr<MeshLib::Mesh> const result(MeshToolsLib::addLayerToMesh(
        *mesh, height, "mesh", false, copy_material_ids, std::nullopt));

    ASSERT_EQ(2 * mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(2 * mesh->getNumberOfElements(), result->getNumberOfElements());

    auto const& n_elems =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*result);
    ASSERT_EQ(50, n_elems.at(MeshLib::MeshElemType::TRIANGLE));
    ASSERT_EQ(50, n_elems.at(MeshLib::MeshElemType::PRISM));

    AddLayerValidation::testZCoords2D(*mesh, *result, -1 * height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddTopLayerToQuadMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshToolsLib::MeshGenerator::generateRegularQuadMesh(5, 5));
    double const height(1);
    constexpr bool copy_material_ids = false;
    std::unique_ptr<MeshLib::Mesh> const result(MeshToolsLib::addLayerToMesh(
        *mesh, height, "mesh", true, copy_material_ids, std::nullopt));

    ASSERT_EQ(2 * mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(2 * mesh->getNumberOfElements(), result->getNumberOfElements());

    auto const& n_elems =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*result);
    ASSERT_EQ(25, n_elems.at(MeshLib::MeshElemType::QUAD));
    ASSERT_EQ(25, n_elems.at(MeshLib::MeshElemType::HEXAHEDRON));

    AddLayerValidation::testZCoords2D(*mesh, *result, height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddBottomLayerToQuadMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshToolsLib::MeshGenerator::generateRegularQuadMesh(5, 5));
    double const height(1);
    constexpr bool copy_material_ids = false;
    std::unique_ptr<MeshLib::Mesh> const result(MeshToolsLib::addLayerToMesh(
        *mesh, height, "mesh", false, copy_material_ids, std::nullopt));

    ASSERT_EQ(2 * mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(2 * mesh->getNumberOfElements(), result->getNumberOfElements());

    auto const& n_elems =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*result);
    ASSERT_EQ(25, n_elems.at(MeshLib::MeshElemType::QUAD));
    ASSERT_EQ(25, n_elems.at(MeshLib::MeshElemType::HEXAHEDRON));

    AddLayerValidation::testZCoords2D(*mesh, *result, -1 * height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddTopLayerToHexMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshToolsLib::MeshGenerator::generateRegularHexMesh(5, 5));
    double const height(1);
    constexpr bool copy_material_ids = false;
    std::unique_ptr<MeshLib::Mesh> const result(MeshToolsLib::addLayerToMesh(
        *mesh, height, "mesh", true, copy_material_ids, std::nullopt));

    ASSERT_EQ(mesh->getNumberOfNodes(), result->getNumberOfNodes() - 36);
    ASSERT_EQ(mesh->getNumberOfElements(), result->getNumberOfElements() - 25);

    auto const& n_elems =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*result);
    ASSERT_EQ(150, n_elems.at(MeshLib::MeshElemType::HEXAHEDRON));

    Eigen::Vector3d const dir({0, 0, -1});
    std::unique_ptr<MeshLib::Mesh> const test_input(
        MeshToolsLib::MeshSurfaceExtraction::getMeshSurface(*mesh, dir, 90));
    std::unique_ptr<MeshLib::Mesh> const test_output(
        MeshToolsLib::MeshSurfaceExtraction::getMeshSurface(*result, dir, 90));
    AddLayerValidation::testZCoords3D(*test_input, *test_output, height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddBottomLayerToHexMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshToolsLib::MeshGenerator::generateRegularHexMesh(5, 5));
    double const height(1);
    constexpr bool copy_material_ids = false;
    std::unique_ptr<MeshLib::Mesh> const result(MeshToolsLib::addLayerToMesh(
        *mesh, height, "mesh", false, copy_material_ids, std::nullopt));

    ASSERT_EQ(mesh->getNumberOfNodes(), result->getNumberOfNodes() - 36);
    ASSERT_EQ(mesh->getNumberOfElements(), result->getNumberOfElements() - 25);

    auto const& n_elems =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*result);
    ASSERT_EQ(150, n_elems.at(MeshLib::MeshElemType::HEXAHEDRON));

    Eigen::Vector3d const dir({0, 0, 1});
    std::unique_ptr<MeshLib::Mesh> const test_input(
        MeshToolsLib::MeshSurfaceExtraction::getMeshSurface(*mesh, dir, 90));
    std::unique_ptr<MeshLib::Mesh> const test_output(
        MeshToolsLib::MeshSurfaceExtraction::getMeshSurface(*result, dir, 90));
    AddLayerValidation::testZCoords3D(*test_input, *test_output, -1 * height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddTopLayerToPrismMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshToolsLib::MeshGenerator::generateRegularTriMesh(5, 5));
    constexpr bool copy_material_ids = false;
    std::unique_ptr<MeshLib::Mesh> const mesh2(MeshToolsLib::addLayerToMesh(
        *mesh, 5, "mesh", true, copy_material_ids, std::nullopt));
    double const height(1);
    std::unique_ptr<MeshLib::Mesh> const result(MeshToolsLib::addLayerToMesh(
        *mesh2, height, "mesh", true, copy_material_ids, std::nullopt));

    ASSERT_EQ(mesh2->getNumberOfNodes() / 2.0 * 3, result->getNumberOfNodes());
    ASSERT_EQ(mesh2->getNumberOfElements() / 2.0 * 3,
              result->getNumberOfElements());

    auto const& n_elems =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*result);
    ASSERT_EQ(50, n_elems.at(MeshLib::MeshElemType::TRIANGLE));
    ASSERT_EQ(100, n_elems.at(MeshLib::MeshElemType::PRISM));

    Eigen::Vector3d const dir({0, 0, -1});
    std::unique_ptr<MeshLib::Mesh> test_input(
        MeshToolsLib::MeshSurfaceExtraction::getMeshSurface(*mesh2, dir, 90));
    std::unique_ptr<MeshLib::Mesh> test_output(
        MeshToolsLib::MeshSurfaceExtraction::getMeshSurface(*result, dir, 90));
    AddLayerValidation::testZCoords3D(*test_input, *test_output, height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddBottomLayerToPrismMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshToolsLib::MeshGenerator::generateRegularTriMesh(5, 5));
    constexpr bool copy_material_ids = false;
    std::unique_ptr<MeshLib::Mesh> const mesh2(MeshToolsLib::addLayerToMesh(
        *mesh, 5, "mesh", true, copy_material_ids, std::nullopt));
    double const height(1);
    std::string const& mat_name("MaterialIDs");
    auto* const mats = mesh2->getProperties().createNewPropertyVector<int>(
        mat_name, MeshLib::MeshItemType::Cell, mesh2->getNumberOfElements(), 1);
    if (mats)
    {
        ranges::fill(mats->begin(), mats->end(), 0);
    }

    std::unique_ptr<MeshLib::Mesh> const result(MeshToolsLib::addLayerToMesh(
        *mesh2, height, "mesh", false, copy_material_ids, std::nullopt));
    ASSERT_EQ(mesh2->getNumberOfNodes() / 2.0 * 3, result->getNumberOfNodes());
    ASSERT_EQ(mesh2->getNumberOfElements() / 2.0 * 3,
              result->getNumberOfElements());

    auto const& n_elems =
        MeshToolsLib::MeshInformation::getNumberOfElementTypes(*result);
    ASSERT_EQ(mesh->getNumberOfElements(),
              n_elems.at(MeshLib::MeshElemType::TRIANGLE));
    ASSERT_EQ(2 * mesh->getNumberOfElements(),
              n_elems.at(MeshLib::MeshElemType::PRISM));
    ASSERT_EQ(1, result->getProperties().size());
    auto const* const new_mats =
        result->getProperties().getPropertyVector<int>(mat_name);
    ASSERT_EQ(result->getNumberOfElements(), new_mats->size());
    ASSERT_EQ(mesh2->getNumberOfElements(),
              std::count(new_mats->cbegin(), new_mats->cend(), 0));
    ASSERT_EQ(mesh->getNumberOfElements(),
              std::count(new_mats->cbegin(), new_mats->cend(), 1));

    Eigen::Vector3d const dir({0, 0, 1});
    std::unique_ptr<MeshLib::Mesh> test_input(
        MeshToolsLib::MeshSurfaceExtraction::getMeshSurface(*mesh2, dir, 90));
    std::unique_ptr<MeshLib::Mesh> test_output(
        MeshToolsLib::MeshSurfaceExtraction::getMeshSurface(*result, dir, 90));
    AddLayerValidation::testZCoords3D(*test_input, *test_output, -1 * height);
    AddLayerValidation::validate(*result, true);
}
