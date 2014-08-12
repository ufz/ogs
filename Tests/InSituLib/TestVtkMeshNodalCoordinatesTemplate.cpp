/**
 * \file
 * \author Lars Bilke
 * \date   2014-02-26
 * \brief  Unit tests for In-Situ mesh nodes
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>

#include "gtest/gtest.h"

#include "Mesh.h"
#include "MeshGenerators/MeshGenerator.h"

#include "VtkMeshNodalCoordinatesTemplate.h"
#include "VtkMappedMesh.h"
#include "VtkMappedMeshSource.h"

class InSituMesh : public ::testing::Test
{
	public:
	InSituMesh()
	 : mesh(nullptr)
	{
		mesh = MeshLib::MeshGenerator::generateRegularQuadMesh(this->length, this->subdivisions);
	}

	~InSituMesh()
	{
	    delete mesh;
	}

	MeshLib::Mesh const* mesh;
	const size_t subdivisions = 99;
	const double length = 10.0;
	const double dx = length / subdivisions;
};

TEST_F(InSituMesh, Construction)
{
	ASSERT_TRUE(mesh != nullptr);
	ASSERT_EQ((subdivisions+1)*(subdivisions+1), mesh->getNNodes());
}


TEST_F(InSituMesh, MappedMesh)
{
	ASSERT_TRUE(mesh != nullptr);

	vtkNew<InSituLib::VtkMappedMesh> vtkMesh;
	vtkMesh->GetImplementation()->SetNodes(mesh->getNodes());
	vtkMesh->GetImplementation()->SetElements(mesh->getElements());

	ASSERT_EQ(subdivisions*subdivisions, vtkMesh->GetNumberOfCells());
	ASSERT_EQ(VTK_QUAD, vtkMesh->GetCellType(0));
	ASSERT_EQ(VTK_QUAD, vtkMesh->GetCellType(vtkMesh->GetNumberOfCells()-1));
	ASSERT_EQ(1, vtkMesh->IsHomogeneous());
	ASSERT_EQ(4, vtkMesh->GetMaxCellSize());


	ASSERT_EQ(0, vtkMesh->GetNumberOfPoints()); // No points are defined
}

TEST(InSituLibNodalCoordinates, Init)
{
	const size_t subdivisions = 99;
	const double length = 10.0;
	const double dx = length / subdivisions;

	MeshLib::Mesh* mesh = MeshLib::MeshGenerator::generateRegularQuadMesh(length, subdivisions);

	vtkNew<InSituLib::VtkMeshNodalCoordinatesTemplate<double> > nodeCoords;
	nodeCoords->SetNodes(mesh->getNodes());
	//nodeCoords->PrintSelf(std::cout, vtkIndent());

	ASSERT_EQ(nodeCoords->GetNumberOfComponents(), 3);
	const size_t numTuples = (subdivisions+1)*(subdivisions+1);
	ASSERT_EQ(nodeCoords->GetNumberOfTuples(), numTuples);

	// First point
	ASSERT_EQ(nodeCoords->GetValueReference(0), 0.0);
	ASSERT_EQ(nodeCoords->GetValueReference(1), 0.0);
	ASSERT_EQ(nodeCoords->GetValueReference(2), 0.0);

	// First second
	ASSERT_NEAR(nodeCoords->GetValueReference(3), dx, std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(nodeCoords->GetValueReference(4), 0.0, std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(nodeCoords->GetValueReference(5), 0.0, std::numeric_limits<double>::epsilon());

	double* coords = nodeCoords->GetTuple(1);
	ASSERT_NEAR(coords[0], dx, std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(coords[1], 0.0, std::numeric_limits<double>::epsilon());
	ASSERT_NEAR(coords[2], 0.0, std::numeric_limits<double>::epsilon());

	delete mesh;
}
