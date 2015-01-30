/**
* \file
* \author Lars Bilke
* \date   2014-02-26
* \brief  Unit tests for In-Situ data arrays
*
* \copyright
* Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/
#include <numeric>

#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>

#include "gtest/gtest.h"

#include "Mesh.h"
#include "MeshGenerators/MeshGenerator.h"

#include "VtkMappedElementDataArrayTemplate.h"

TEST(InSituLibMappedArrays, Double)
{
	const size_t mesh_size = 5;
	const double length = 1.0;

	MeshLib::Mesh* mesh = MeshLib::MeshGenerator::generateRegularHexMesh(length, mesh_size);

	ASSERT_TRUE(mesh != nullptr);
	const std::size_t size(mesh_size*mesh_size*mesh_size);

	std::string const prop_name("TestProperty");
	boost::optional<MeshLib::PropertyVector<double> &> double_properties(
		mesh->getProperties().createNewPropertyVector<double>(prop_name,
			MeshLib::MeshItemType::Cell));
	(*double_properties).resize(size);
	std::iota((*double_properties).begin(), (*double_properties).end(), 1);

	vtkNew<InSituLib::VtkMappedElementDataArrayTemplate<double> > dataArray;
	dataArray->SetPropertyVector(*double_properties);

	ASSERT_EQ(dataArray->GetNumberOfComponents(), 1);
	ASSERT_EQ(dataArray->GetNumberOfTuples(), size);

	ASSERT_EQ(dataArray->GetValueReference(0), 1.0);
	double* range = dataArray->GetRange(0);
	ASSERT_EQ(range[0], 1.0);
	ASSERT_EQ(range[1], 1.0 + mesh->getNElements() - 1.0);

	delete mesh;
}

TEST(InSituLibMappedArrays, Int)
{
	const size_t mesh_size = 5;
	const double length = 1.0;

	MeshLib::Mesh* mesh = MeshLib::MeshGenerator::generateRegularHexMesh(length, mesh_size);

	ASSERT_TRUE(mesh != nullptr);
	const std::size_t size(mesh_size*mesh_size*mesh_size);

	std::string const prop_name("TestProperty");
	boost::optional<MeshLib::PropertyVector<int> &> properties(
		mesh->getProperties().createNewPropertyVector<int>(prop_name,
			MeshLib::MeshItemType::Cell));
	(*properties).resize(size);
	std::iota((*properties).begin(), (*properties).end(), 1);

	vtkNew<InSituLib::VtkMappedElementDataArrayTemplate<int> > dataArray;
	dataArray->SetPropertyVector(*properties);

	ASSERT_EQ(dataArray->GetNumberOfComponents(), 1);
	ASSERT_EQ(dataArray->GetNumberOfTuples(), size);

	ASSERT_EQ(dataArray->GetValueReference(0), 1);
	double* range = dataArray->GetRange(0);
	ASSERT_EQ(range[0], 1);
	ASSERT_EQ(range[1], 1 + mesh->getNElements() - 1);

	delete mesh;
}