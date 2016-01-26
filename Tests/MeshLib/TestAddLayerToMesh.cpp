/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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
		int reduce_tests = (testNodeOrder) ? 0 : 1;

		const std::size_t nErrorFlags (static_cast<std::size_t>(ElementErrorFlag::MaxValue));
		ElementErrorFlag flags[nErrorFlags] = {ElementErrorFlag::ZeroVolume, 
		ElementErrorFlag::NonCoplanar, ElementErrorFlag::NonConvex,  ElementErrorFlag::NodeOrder};
		std::vector<ElementErrorCode> codes (MeshLib::MeshValidation::testElementGeometry(mesh));
		for (std::size_t i=0; i<codes.size(); ++i)
			for (std::size_t j=0; j<nErrorFlags-reduce_tests; ++j)
				ASSERT_EQ(false, codes[i][flags[j]]);
	}

	void testZCoords2D(MeshLib::Mesh const& input, MeshLib::Mesh const& output, int height)
	{
		std::size_t nNodes (input.getNNodes());
		for (std::size_t i=0; i<nNodes; ++i)
		{
			ASSERT_EQ((*input.getNode(i))[2], (*output.getNode(i))[2]);
			ASSERT_EQ((*input.getNode(i))[2] + height, (*output.getNode(nNodes+i))[2]);
		}
	}

	void testZCoords3D(MeshLib::Mesh const& input, MeshLib::Mesh const& output, int height)
	{
		std::size_t nNodes (input.getNNodes());
		for (std::size_t i=0; i<nNodes; ++i)
			ASSERT_EQ((*input.getNode(i))[2] + height, (*output.getNode(i))[2]);
	}
};

TEST(MeshLib, AddLayerToLineMesh)
{
	std::unique_ptr<MeshLib::Mesh> mesh (MeshLib::MeshGenerator::generateLineMesh(1.0, 5));

	// add top layer
	int height (1);
	std::unique_ptr<MeshLib::Mesh> result (MeshLib::addLayerToMesh(*mesh, height, "mesh", true));

	ASSERT_EQ(2*mesh->getNNodes(), result->getNNodes());
	ASSERT_EQ(2*mesh->getNElements(), result->getNElements());

	std::array<unsigned, 7> n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
	ASSERT_EQ(5, n_elems[0]);
	ASSERT_EQ(5, n_elems[2]);

	AddLayerValidation::testZCoords2D(*mesh, *result, height);
	AddLayerValidation::validate(*result, false);

	// add bottom layer
	result.reset(MeshLib::addLayerToMesh(*mesh, height, "mesh", false));

	ASSERT_EQ(2*mesh->getNNodes(), result->getNNodes());
	ASSERT_EQ(2*mesh->getNElements(), result->getNElements());

	n_elems = MeshLib::MeshInformation::getNumberOfElementTypes(*result);
	ASSERT_EQ(5, n_elems[0]);
	ASSERT_EQ(5, n_elems[2]);

	AddLayerValidation::testZCoords2D(*mesh, *result, -1 * height);
	AddLayerValidation::validate(*result, false);
}

TEST(MeshLib, AddLayerToTriMesh)
{
	std::unique_ptr<MeshLib::Mesh> mesh (MeshLib::MeshGenerator::generateRegularTriMesh(5, 5));

	// add top layer
	int height (1);
	std::unique_ptr<MeshLib::Mesh> result (MeshLib::addLayerToMesh(*mesh, height, "mesh", true));

	ASSERT_EQ(2*mesh->getNNodes(), result->getNNodes());
	ASSERT_EQ(2*mesh->getNElements(), result->getNElements());

	std::array<unsigned, 7> n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
	ASSERT_EQ(50, n_elems[1]);
	ASSERT_EQ(50, n_elems[6]);

	AddLayerValidation::testZCoords2D(*mesh, *result, height);
	AddLayerValidation::validate(*result, true);

	// add bottom layer
	result.reset(MeshLib::addLayerToMesh(*mesh, height, "mesh", false));

	ASSERT_EQ(2*mesh->getNNodes(), result->getNNodes());
	ASSERT_EQ(2*mesh->getNElements(), result->getNElements());

	n_elems = MeshLib::MeshInformation::getNumberOfElementTypes(*result);
	ASSERT_EQ(50, n_elems[1]);
	ASSERT_EQ(50, n_elems[6]);

	AddLayerValidation::testZCoords2D(*mesh, *result, -1 * height);
	AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddLayerToQuadMesh)
{
	std::unique_ptr<MeshLib::Mesh> mesh (MeshLib::MeshGenerator::generateRegularQuadMesh(5, 5));

	// add top layer
	int height (1);
	std::unique_ptr<MeshLib::Mesh> result (MeshLib::addLayerToMesh(*mesh, height, "mesh", true));

	ASSERT_EQ(2*mesh->getNNodes(), result->getNNodes());
	ASSERT_EQ(2*mesh->getNElements(), result->getNElements());

	std::array<unsigned, 7> n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
	ASSERT_EQ(25, n_elems[2]);
	ASSERT_EQ(25, n_elems[4]);

	AddLayerValidation::testZCoords2D(*mesh, *result, height);
	AddLayerValidation::validate(*result, true);

	// add bottom layer
	result.reset(MeshLib::addLayerToMesh(*mesh, height, "mesh", false));

	ASSERT_EQ(2*mesh->getNNodes(), result->getNNodes());
	ASSERT_EQ(2*mesh->getNElements(), result->getNElements());

	n_elems = MeshLib::MeshInformation::getNumberOfElementTypes(*result);
	ASSERT_EQ(25, n_elems[2]);
	ASSERT_EQ(25, n_elems[4]);

	AddLayerValidation::testZCoords2D(*mesh, *result, -1 * height);
	AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddLayerToHexMesh)
{
	std::unique_ptr<MeshLib::Mesh> mesh (MeshLib::MeshGenerator::generateRegularHexMesh(5, 5));

	// add top layer
	int height (1);
	std::unique_ptr<MeshLib::Mesh> result (MeshLib::addLayerToMesh(*mesh, height, "mesh", true));

	ASSERT_EQ(mesh->getNNodes(), result->getNNodes()-36);
	ASSERT_EQ(mesh->getNElements(), result->getNElements()-25);

	std::array<unsigned, 7> n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
	ASSERT_EQ(150, n_elems[4]);

	MathLib::Vector3 const dir(0, 0, -1);
	std::unique_ptr<MeshLib::Mesh> test_input (
		MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh, dir, 90));
	std::unique_ptr<MeshLib::Mesh> test_output (
		MeshLib::MeshSurfaceExtraction::getMeshSurface(*result, dir, 90));
	AddLayerValidation::testZCoords3D(*test_input, *test_output, height);
	AddLayerValidation::validate(*result, true);

	// add bottom layer
	result.reset(MeshLib::addLayerToMesh(*mesh, height, "mesh", false));

	ASSERT_EQ(mesh->getNNodes(), result->getNNodes()-36);
	ASSERT_EQ(mesh->getNElements(), result->getNElements()-25);

	n_elems = MeshLib::MeshInformation::getNumberOfElementTypes(*result);
	ASSERT_EQ(150, n_elems[4]);

	MathLib::Vector3 const dir2(0, 0, 1);
	test_input.reset (MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh, dir2, 90));
	test_output.reset (MeshLib::MeshSurfaceExtraction::getMeshSurface(*result, dir2, 90));
	AddLayerValidation::testZCoords3D(*test_input, *test_output, -1 * height);
	AddLayerValidation::validate(*result, true);
}

TEST(MeshLib, AddLayerToPrismMesh)
{
	std::unique_ptr<MeshLib::Mesh> mesh (MeshLib::MeshGenerator::generateRegularTriMesh(5, 5));
	// create prism-mesh from tri-mesh
	std::unique_ptr<MeshLib::Mesh> mesh2 (MeshLib::addLayerToMesh(*mesh, 5, "mesh", true));

	// add top layer
	int height (1);
	std::unique_ptr<MeshLib::Mesh> result (MeshLib::addLayerToMesh(*mesh2, height, "mesh", true));

	ASSERT_EQ(mesh2->getNNodes()/2.0 * 3, result->getNNodes());
	ASSERT_EQ(mesh2->getNElements()/2.0 * 3, result->getNElements());

	std::array<unsigned, 7> n_elems (MeshLib::MeshInformation::getNumberOfElementTypes(*result));
	ASSERT_EQ(50, n_elems[1]);
	ASSERT_EQ(100, n_elems[6]);

	MathLib::Vector3 const dir(0, 0, -1);
	std::unique_ptr<MeshLib::Mesh> test_input (
		MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh2, dir, 90));
	std::unique_ptr<MeshLib::Mesh> test_output (
		MeshLib::MeshSurfaceExtraction::getMeshSurface(*result, dir, 90));
	AddLayerValidation::testZCoords3D(*test_input, *test_output, height);
	AddLayerValidation::validate(*result, true);

	// add bottom layer
	result.reset(MeshLib::addLayerToMesh(*mesh2, height, "mesh", false));

	ASSERT_EQ(mesh2->getNNodes()/2.0 * 3, result->getNNodes());
	ASSERT_EQ(mesh2->getNElements()/2.0 * 3, result->getNElements());

	n_elems = MeshLib::MeshInformation::getNumberOfElementTypes(*result);
	ASSERT_EQ(50, n_elems[1]);
	ASSERT_EQ(100, n_elems[6]);

	MathLib::Vector3 const dir2(0, 0, 1);
	test_input.reset (MeshLib::MeshSurfaceExtraction::getMeshSurface(*mesh2, dir2, 90));
	test_output.reset (MeshLib::MeshSurfaceExtraction::getMeshSurface(*result, dir2, 90));
	AddLayerValidation::testZCoords3D(*test_input, *test_output, -1 * height);
	AddLayerValidation::validate(*result, true);
}

