/**
* \file
* \author Lars Bilke
* \date   2014-02-26
* \brief  Unit tests for In-Situ data arrays
*
* \copyright
* Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
*            Distributed under a Modified BSD License.
*              See accompanying file LICENSE.txt or
*              http://www.opengeosys.org/project/license
*
*/
#include <numeric>

#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnsignedIntArray.h>

#include "gtest/gtest.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/Vtk/VtkMappedPropertyVectorTemplate.h"

// Creates a PropertyVector<double> and maps it into a vtkDataArray-equivalent
TEST(MeshLibMappedPropertyVector, Double)
{
    const std::size_t mesh_size = 5;
    const double length = 1.0;

    MeshLib::Mesh* mesh = MeshLib::MeshGenerator::generateRegularHexMesh(length, mesh_size);

    ASSERT_TRUE(mesh != nullptr);
    const std::size_t number_of_tuples(mesh_size*mesh_size*mesh_size);

    std::string const prop_name("TestProperty");
    auto* const double_properties =
        mesh->getProperties().createNewPropertyVector<double>(
            prop_name, MeshLib::MeshItemType::Cell);
    double_properties->resize(number_of_tuples);
    std::iota(double_properties->begin(), double_properties->end(), 1);

    vtkNew<vtkDoubleArray> dataArray;
    dataArray->SetNumberOfComponents(1);
    dataArray->SetArray(double_properties->data(),
        static_cast<vtkIdType>(double_properties->size()), 1);

    ASSERT_EQ(dataArray->GetNumberOfComponents(), 1);
    ASSERT_EQ(dataArray->GetNumberOfTuples(), number_of_tuples);

    ASSERT_EQ(dataArray->GetValue(0), 1.0);
    double* range = dataArray->GetRange(0);
    ASSERT_EQ(range[0], 1.0);
    ASSERT_EQ(range[1], 1.0 + mesh->getNumberOfElements() - 1.0);

    delete mesh;
}

// Creates a PropertyVector<int> and maps it into a vtkDataArray-equivalent
TEST(MeshLibMappedPropertyVector, Int)
{
    const std::size_t mesh_size = 5;
    const double length = 1.0;

    MeshLib::Mesh* mesh = MeshLib::MeshGenerator::generateRegularHexMesh(length, mesh_size);

    ASSERT_TRUE(mesh != nullptr);
    const std::size_t number_of_tuples(mesh_size*mesh_size*mesh_size);

    std::string const prop_name("TestProperty");
    auto* const properties = mesh->getProperties().createNewPropertyVector<int>(
        prop_name, MeshLib::MeshItemType::Cell);
    properties->resize(number_of_tuples);
    std::iota(properties->begin(), properties->end(), -5);

    vtkNew<vtkIntArray> dataArray;
    dataArray->SetNumberOfComponents(1);
    dataArray->SetArray(properties->data(),
        static_cast<vtkIdType>(properties->size()), 1);

    ASSERT_EQ(dataArray->GetNumberOfComponents(), 1);
    ASSERT_EQ(dataArray->GetNumberOfTuples(), number_of_tuples);

    ASSERT_EQ(dataArray->GetValue(0), -5);
    double* range = dataArray->GetRange(0);
    ASSERT_EQ(-5.0, range[0]);
    ASSERT_EQ(-5.0 + static_cast<double>(mesh->getNumberOfElements()) - 1.0, range[1]);

    delete mesh;
}

// Creates a PropertyVector<unsigned> and maps it into a vtkDataArray-equivalent
TEST(MeshLibMappedPropertyVector, Unsigned)
{
    const std::size_t mesh_size = 5;
    const double length = 1.0;

    MeshLib::Mesh* mesh = MeshLib::MeshGenerator::generateRegularHexMesh(length, mesh_size);

    ASSERT_TRUE(mesh != nullptr);
    const std::size_t number_of_tuples(mesh_size*mesh_size*mesh_size);

    std::string const prop_name("TestProperty");
    auto* const properties =
        mesh->getProperties().createNewPropertyVector<unsigned>(
            prop_name, MeshLib::MeshItemType::Cell);
    properties->resize(number_of_tuples);
    std::iota(properties->begin(), properties->end(), 0);

    vtkNew<vtkUnsignedIntArray> dataArray;
    dataArray->SetNumberOfComponents(1);
    dataArray->SetArray(properties->data(),
        static_cast<vtkIdType>(properties->size()), 1);

    ASSERT_EQ(dataArray->GetNumberOfComponents(), 1);
    ASSERT_EQ(dataArray->GetNumberOfTuples(), number_of_tuples);

    ASSERT_EQ(dataArray->GetValue(0), 0);
    double* range = dataArray->GetRange(0);
    ASSERT_EQ(range[0], 0);
    ASSERT_EQ(range[1], 0 + mesh->getNumberOfElements() - 1);

    delete mesh;
}
