/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>

#include "gtest/gtest.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/MeshSurfaceExtraction.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshEditing/AddLayerToMesh.h"
#include "MeshLib/MeshQuality/MeshValidation.h"

namespace AddLayerValidation
{
    // validates mesh. for line meshes, node order tests will fail because vertical elements are
    // created during the extrusion. Therefore node order tests are switched off for line meshes.
    void validate(MeshLib::Mesh const& mesh, bool testNodeOrder)
    {
        int const reduce_tests = (testNodeOrder) ? 0 : 1;

        auto const nErrorFlags(
            static_cast<std::size_t>(ElementErrorFlag::MaxValue));
        ElementErrorFlag const flags[nErrorFlags] = {ElementErrorFlag::ZeroVolume,
        ElementErrorFlag::NonCoplanar, ElementErrorFlag::NonConvex,  ElementErrorFlag::NodeOrder};
        std::vector<ElementErrorCode> const codes (MeshLib::MeshValidation::testElementGeometry(mesh));
        for (auto code : codes)
            for (std::size_t j=0; j<nErrorFlags-reduce_tests; ++j)
                ASSERT_FALSE(code[flags[j]]);
    }

    void testZCoords2D(MeshLib::Mesh const& input, MeshLib::Mesh const& output, double height)
    {
        std::size_t const nNodes (input.getNumberOfNodes());
        for (std::size_t i=0; i<nNodes; ++i)
        {
            ASSERT_EQ((*input.getNode(i))[2], (*output.getNode(i))[2]);
            ASSERT_EQ((*input.getNode(i))[2] + height, (*output.getNode(nNodes+i))[2]);
        }
    }

    void testZCoords3D(MeshLib::Mesh const& input, MeshLib::Mesh const& output, double height)
    {
        std::size_t const nNodes (input.getNumberOfNodes());
        for (std::size_t i=0; i<nNodes; ++i)
            ASSERT_EQ((*input.getNode(i))[2] + height, (*output.getNode(i))[2]);
    }
};

TEST(MeshLib, AddTopLayerToLineMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh (MeshLib::MeshGenerator::generateLineMesh(1.0, 5));
    double const height (1);
    std::unique_ptr<MeshLib::Mesh> const result (MeshLib::addLayerToMesh(*mesh, height, "mesh", true));

    ASSERT_EQ(2*mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(2*mesh->getNumberOfElements(), result->getNumberOfElements());

    std::array<unsigned, 7> const n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
    ASSERT_EQ(5, n_elems[0]); // tests if 5 lines are present
    ASSERT_EQ(5, n_elems[2]); // tests if 5 quads are present

    AddLayerValidation::testZCoords2D(*mesh, *result, height);
    AddLayerValidation::validate(*result, false);
}

TEST(MeshLib, AddBottomLayerToLineMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh (MeshLib::MeshGenerator::generateLineMesh(1.0, 5));
    double const height (1);
    std::unique_ptr<MeshLib::Mesh> const result (MeshLib::addLayerToMesh(*mesh, height, "mesh", false));

    ASSERT_EQ(2*mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(2*mesh->getNumberOfElements(), result->getNumberOfElements());

    std::array<unsigned, 7> const n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
    ASSERT_EQ(5, n_elems[0]); // tests if 5 lines are present
    ASSERT_EQ(5, n_elems[2]); // tests if 5 quads are present

    AddLayerValidation::testZCoords2D(*mesh, *result, -1 * height);
    AddLayerValidation::validate(*result, false);
}

TEST(MeshLib, AddTopLayerToTriMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh (MeshLib::MeshGenerator::generateRegularTriMesh(5, 5));
    std::string const& mat_name ("MaterialIDs");
    auto* const mats = mesh->getProperties().createNewPropertyVector<int>(
        mat_name, MeshLib::MeshItemType::Cell);
    if (mats)
    {
        mats->resize(mesh->getNumberOfElements(), 0);
    }
    auto* const test = mesh->getProperties().createNewPropertyVector<double>(
        "test", MeshLib::MeshItemType::Cell);
    if (test)
    {
        test->resize(mesh->getNumberOfElements(), 0.1);
    }
    ASSERT_EQ(2, mesh->getProperties().getPropertyVectorNames().size());
    double const height (1);
    std::unique_ptr<MeshLib::Mesh> const result (MeshLib::addLayerToMesh(*mesh, height, "mesh", true));

    ASSERT_EQ(2*mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(2*mesh->getNumberOfElements(), result->getNumberOfElements());

    std::array<unsigned, 7> const n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
    ASSERT_EQ(mesh->getNumberOfElements(), n_elems[1]); // tests if 50 tris are present
    ASSERT_EQ(mesh->getNumberOfElements(), n_elems[6]); // tests if 50 prisms are present

    ASSERT_EQ(1, result->getProperties().getPropertyVectorNames().size());
    auto const* const new_mats =
        result->getProperties().getPropertyVector<int>(mat_name);
    ASSERT_EQ(result->getNumberOfElements(), new_mats->size());
    ASSERT_EQ(mesh->getNumberOfElements(), std::count(new_mats->cbegin(), new_mats->cend(), 0));
    ASSERT_EQ(mesh->getNumberOfElements(), std::count(new_mats->cbegin(), new_mats->cend(), 1));
    AddLayerValidation::testZCoords2D(*mesh, *result, height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddBottomLayerToTriMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh (MeshLib::MeshGenerator::generateRegularTriMesh(5, 5));
    double const height (1);
    std::unique_ptr<MeshLib::Mesh> const result (MeshLib::addLayerToMesh(*mesh, height, "mesh", false));

    ASSERT_EQ(2*mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(2*mesh->getNumberOfElements(), result->getNumberOfElements());

    std::array<unsigned, 7> const n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
    ASSERT_EQ(50, n_elems[1]); // tests if 50 tris are present
    ASSERT_EQ(50, n_elems[6]); // tests if 50 prisms are present

    AddLayerValidation::testZCoords2D(*mesh, *result, -1 * height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddTopLayerToQuadMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh (MeshLib::MeshGenerator::generateRegularQuadMesh(5, 5));
    double const height (1);
    std::unique_ptr<MeshLib::Mesh> const result (MeshLib::addLayerToMesh(*mesh, height, "mesh", true));

    ASSERT_EQ(2*mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(2*mesh->getNumberOfElements(), result->getNumberOfElements());

    std::array<unsigned, 7> const n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
    ASSERT_EQ(25, n_elems[2]); // tests if 25 quads are present
    ASSERT_EQ(25, n_elems[4]); // tests if 25 hexes are present

    AddLayerValidation::testZCoords2D(*mesh, *result, height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddBottomLayerToQuadMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh (MeshLib::MeshGenerator::generateRegularQuadMesh(5, 5));
    double const height (1);
    std::unique_ptr<MeshLib::Mesh> const result (MeshLib::addLayerToMesh(*mesh, height, "mesh", false));

    ASSERT_EQ(2*mesh->getNumberOfNodes(), result->getNumberOfNodes());
    ASSERT_EQ(2*mesh->getNumberOfElements(), result->getNumberOfElements());

    std::array<unsigned, 7> const n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
    ASSERT_EQ(25, n_elems[2]); // tests if 25 quads are present
    ASSERT_EQ(25, n_elems[4]); // tests if 25 hexes are present

    AddLayerValidation::testZCoords2D(*mesh, *result, -1 * height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddTopLayerToHexMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh (MeshLib::MeshGenerator::generateRegularHexMesh(5, 5));
    double const height (1);
    std::unique_ptr<MeshLib::Mesh> const result (MeshLib::addLayerToMesh(*mesh, height, "mesh", true));

    ASSERT_EQ(mesh->getNumberOfNodes(), result->getNumberOfNodes()-36);
    ASSERT_EQ(mesh->getNumberOfElements(), result->getNumberOfElements()-25);

    std::array<unsigned, 7> const n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
    ASSERT_EQ(150, n_elems[4]); // tests if 150 hexes are present

    MathLib::Vector3 const dir(0, 0, -1);
    std::unique_ptr<MeshLib::Mesh> const test_input (
        MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh, dir, 90));
    std::unique_ptr<MeshLib::Mesh> const test_output (
        MeshLib::MeshSurfaceExtraction::getMeshSurface(*result, dir, 90));
    AddLayerValidation::testZCoords3D(*test_input, *test_output, height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddBottomLayerToHexMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh (MeshLib::MeshGenerator::generateRegularHexMesh(5, 5));
    double const height (1);
    std::unique_ptr<MeshLib::Mesh> const result (MeshLib::addLayerToMesh(*mesh, height, "mesh", false));

    ASSERT_EQ(mesh->getNumberOfNodes(), result->getNumberOfNodes()-36);
    ASSERT_EQ(mesh->getNumberOfElements(), result->getNumberOfElements()-25);

    std::array<unsigned, 7> const n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
    ASSERT_EQ(150, n_elems[4]); // tests if 150 hexes are present

    MathLib::Vector3 const dir(0, 0, 1);
    std::unique_ptr<MeshLib::Mesh> const test_input (
        MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh, dir, 90));
    std::unique_ptr<MeshLib::Mesh> const test_output (
        MeshLib::MeshSurfaceExtraction::getMeshSurface(*result, dir, 90));
    AddLayerValidation::testZCoords3D(*test_input, *test_output, -1 * height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddTopLayerToPrismMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh (MeshLib::MeshGenerator::generateRegularTriMesh(5, 5));
    std::unique_ptr<MeshLib::Mesh> const mesh2 (MeshLib::addLayerToMesh(*mesh, 5, "mesh", true));
    double const height (1);
    std::unique_ptr<MeshLib::Mesh> const result (MeshLib::addLayerToMesh(*mesh2, height, "mesh", true));

    ASSERT_EQ(mesh2->getNumberOfNodes()/2.0 * 3, result->getNumberOfNodes());
    ASSERT_EQ(mesh2->getNumberOfElements()/2.0 * 3, result->getNumberOfElements());

    std::array<unsigned, 7> const n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
    ASSERT_EQ(50, n_elems[1]); // tests if 50 tris are present
    ASSERT_EQ(100, n_elems[6]); // tests if 50 prisms are present

    MathLib::Vector3 const dir(0, 0, -1);
    std::unique_ptr<MeshLib::Mesh> test_input (
        MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh2, dir, 90));
    std::unique_ptr<MeshLib::Mesh> test_output (
        MeshLib::MeshSurfaceExtraction::getMeshSurface(*result, dir, 90));
    AddLayerValidation::testZCoords3D(*test_input, *test_output, height);
    AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddBottomLayerToPrismMesh)
{
    std::unique_ptr<MeshLib::Mesh> const mesh (MeshLib::MeshGenerator::generateRegularTriMesh(5, 5));
    std::unique_ptr<MeshLib::Mesh> const mesh2 (MeshLib::addLayerToMesh(*mesh, 5, "mesh", true));
    double const height (1);
    std::string const& mat_name ("MaterialIDs");
    auto* const mats = mesh2->getProperties().createNewPropertyVector<int>(
        mat_name, MeshLib::MeshItemType::Cell);
    if (mats)
    {
        mats->resize(mesh2->getNumberOfElements(), 0);
    }

    std::unique_ptr<MeshLib::Mesh> const result (MeshLib::addLayerToMesh(*mesh2, height, "mesh", false));
    ASSERT_EQ(mesh2->getNumberOfNodes()/2.0 * 3, result->getNumberOfNodes());
    ASSERT_EQ(mesh2->getNumberOfElements()/2.0 * 3, result->getNumberOfElements());

    std::array<unsigned, 7> const n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
    ASSERT_EQ(mesh->getNumberOfElements(), n_elems[1]); // tests if 50 tris are present
    ASSERT_EQ(2 * mesh->getNumberOfElements(), n_elems[6]); // tests if 50 prisms are present
    ASSERT_EQ(1, result->getProperties().getPropertyVectorNames().size());
    auto const* const new_mats =
        result->getProperties().getPropertyVector<int>(mat_name);
    ASSERT_EQ(result->getNumberOfElements(), new_mats->size());
    ASSERT_EQ(mesh2->getNumberOfElements(), std::count(new_mats->cbegin(), new_mats->cend(), 0));
    ASSERT_EQ(mesh->getNumberOfElements(), std::count(new_mats->cbegin(), new_mats->cend(), 1));

    MathLib::Vector3 const dir(0, 0, 1);
    std::unique_ptr<MeshLib::Mesh> test_input (
        MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh2, dir, 90));
    std::unique_ptr<MeshLib::Mesh> test_output (
        MeshLib::MeshSurfaceExtraction::getMeshSurface(*result, dir, 90));
    AddLayerValidation::testZCoords3D(*test_input, *test_output, -1 * height);
    AddLayerValidation::validate(*result, true);
}
