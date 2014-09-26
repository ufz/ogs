/**
 * \file
 * \author Lars Bilke
 * \date   2014-08-12
 * \brief  Unit tests for In-Situ mesh source
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "gtest/gtest.h"

#include "BaseLib/BuildInfo.h"
#include "Mesh.h"
#include "MeshGenerators/MeshGenerator.h"
#include "MeshGenerators/VtkMeshConverter.h"

#include "VtkMappedMesh.h"
#include "VtkMappedMeshSource.h"

#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>

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

TEST_F(InSituMesh, MappedMeshSourceRoundtrip)
{
	// TODO Add more comparison criteria

	ASSERT_TRUE(mesh != nullptr);
	std::string test_data_file(BaseLib::BuildInfo::tests_tmp_path + "/MappedMeshSourceRoundtrip.vtu");

	// Test VtkMappedMeshSource, i.e. OGS mesh to VTK mesh
	vtkNew<InSituLib::VtkMappedMeshSource> vtkSource;
	vtkSource->SetMesh(mesh);
	vtkSource->Update();
	vtkUnstructuredGrid* output = vtkSource->GetOutput();

	ASSERT_EQ((subdivisions+1)*(subdivisions+1), output->GetNumberOfPoints());
	ASSERT_EQ(subdivisions*subdivisions, output->GetNumberOfCells());

	// Write VTK mesh to file
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtuWriter =
		vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	// Setting binary file  mode, otherwise corrupted output due to VTK bug
	// See http://www.paraview.org/Bug/view.php?id=13382
	vtuWriter->SetDataModeToBinary();
	vtuWriter->SetFileName(test_data_file.c_str());
	vtuWriter->SetInputConnection(vtkSource->GetOutputPort());
	vtuWriter->Write();

	// Read back VTK mesh
	vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
		vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
	reader->SetFileName(test_data_file.c_str());
	reader->Update();
	vtkUnstructuredGrid* vtkMesh = reader->GetOutput();

	// Both VTK meshes should be identical
	ASSERT_EQ(vtkMesh->GetNumberOfPoints(), output->GetNumberOfPoints());
	ASSERT_EQ(vtkMesh->GetNumberOfCells(), output->GetNumberOfCells());

	// Both OGS meshes should be identical
	MeshLib::Mesh* newMesh = MeshLib::VtkMeshConverter::convertUnstructuredGrid(vtkMesh);
	ASSERT_EQ(mesh->getNNodes(), newMesh->getNNodes());
	ASSERT_EQ(mesh->getNElements(), newMesh->getNElements());
}
