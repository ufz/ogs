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

#include "BaseLib/BuildInfo.h"

#include "FileIO/VtkIO/VtuInterface.h"

#include "InSituLib/VtkMappedMesh.h"
#include "InSituLib/VtkMappedMeshSource.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshGenerators/VtkMeshConverter.h"

#include "gtest/gtest.h"

#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkCellData.h>

class InSituMesh : public ::testing::Test
{
	public:
	InSituMesh()
	 : mesh(nullptr)
	{
		mesh = MeshLib::MeshGenerator::generateRegularQuadMesh(this->length, this->subdivisions);
		for(MeshLib::Element* element: mesh->getElements())
			element->setValue(0);
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

	// -- Test VtkMappedMeshSource, i.e. OGS mesh to VTK mesh
	vtkNew<InSituLib::VtkMappedMeshSource> vtkSource;
	vtkSource->SetMesh(mesh);
	vtkSource->Update();
	vtkUnstructuredGrid* output = vtkSource->GetOutput();

	// Point and cell numbers
	ASSERT_EQ((subdivisions+1)*(subdivisions+1), output->GetNumberOfPoints());
	ASSERT_EQ(subdivisions*subdivisions, output->GetNumberOfCells());

	// Cell data array
	vtkDataArray* matIdsArray = output->GetCellData()->GetScalars("MaterialIDs");
	ASSERT_EQ(matIdsArray->GetSize(), mesh->getNElements());
	ASSERT_EQ((unsigned)matIdsArray->GetComponent(0, 0), 0);
	double* range = matIdsArray->GetRange(0);
	ASSERT_EQ((unsigned)range[0], 0);
	ASSERT_EQ((unsigned)range[1], 0);

	// -- Write VTK mesh to file (in all combinations of binary, appended and compressed)
	// atm vtkXMLWriter::Appended does not work, see http://www.paraview.org/Bug/view.php?id=13382
	std::array<int, 2> dataModes = {{ vtkXMLWriter::Ascii, vtkXMLWriter::Binary }};
	std::array<bool, 2> booleans = {{ true, false }};
	for(int dataMode : dataModes)
	{
		for(bool compressed : booleans)
		{
			FileIO::VtuInterface vtuInterface(mesh, dataMode, compressed);
			ASSERT_EQ(vtuInterface.writeToFile(test_data_file), 1);

			// -- Read back VTK mesh
			vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
				vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
			reader->SetFileName(test_data_file.c_str());
			reader->Update();
			vtkUnstructuredGrid* vtkMesh = reader->GetOutput();

			// Both VTK meshes should be identical
			ASSERT_EQ(vtkMesh->GetNumberOfPoints(), output->GetNumberOfPoints());
			ASSERT_EQ(vtkMesh->GetNumberOfCells(), output->GetNumberOfCells());
			ASSERT_EQ(vtkMesh->GetCellData()->GetScalars("MaterialIDs")->GetNumberOfTuples(), matIdsArray->GetNumberOfTuples());

			// Both OGS meshes should be identical
			MeshLib::Mesh* newMesh = MeshLib::VtkMeshConverter::convertUnstructuredGrid(vtkMesh);
			ASSERT_EQ(mesh->getNNodes(), newMesh->getNNodes());
			ASSERT_EQ(mesh->getNElements(), newMesh->getNElements());
		}
	}
}
