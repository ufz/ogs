/**
 * \file
 * \author Lars Bilke
 * \date   2014-08-12
 * \brief  Unit tests for In-Situ mesh source
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <gtest/gtest.h>
#include <vtkCellData.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <memory>
#include <numeric>

#include "BaseLib/StringTools.h"
#include "InfoLib/TestInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/VtkIO/VtkMeshConverter.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Vtk/VtkMappedMeshSource.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"

// Creates a mesh with different types of data (double, int, etc.) and point,
// cell, or integration point properties.
class InSituMesh : public ::testing::Test
{
public:
    InSituMesh()
    {
        mesh = MeshToolsLib::MeshGenerator::generateRegularHexMesh(
            this->length, this->subdivisions);

        std::string const point_prop_name("PointDoubleProperty");
        auto* const point_double_properties =
            mesh->getProperties().createNewPropertyVector<double>(
                point_prop_name, MeshLib::MeshItemType::Node,
                mesh->getNumberOfNodes(), 1);
        std::iota(point_double_properties->begin(),
                  point_double_properties->end(),
                  1);

        std::string const cell_prop_name("CellDoubleProperty");
        auto* const cell_double_properties =
            mesh->getProperties().createNewPropertyVector<double>(
                cell_prop_name, MeshLib::MeshItemType::Cell,
                mesh->getNumberOfElements(), 1);
        std::iota(
            cell_double_properties->begin(), cell_double_properties->end(), 1);

        std::string const field_prop_name("FieldDoubleProperty");
        auto* const field_double_properties =
            mesh->getProperties().createNewPropertyVector<double>(
                field_prop_name, MeshLib::MeshItemType::IntegrationPoint);
        field_double_properties->resize(mesh->getNumberOfElements() * 2);
        std::iota(field_double_properties->begin(),
                  field_double_properties->end(),
                  1);

        std::string const point_int_prop_name("PointIntProperty");
        auto* const point_int_properties =
            mesh->getProperties().createNewPropertyVector<int>(
                point_int_prop_name,
                MeshLib::MeshItemType::Node,
                mesh->getNumberOfNodes(),
                1);
        std::iota(
            point_int_properties->begin(), point_int_properties->end(), 1);

        std::string const cell_int_prop_name("CellIntProperty");
        auto* const cell_int_properties =
            mesh->getProperties().createNewPropertyVector<int>(
                cell_int_prop_name, MeshLib::MeshItemType::Cell,
                mesh->getNumberOfElements(), 1);
        std::iota(cell_int_properties->begin(), cell_int_properties->end(), 1);

        std::string const field_int_prop_name("FieldIntProperty");
        auto* const field_int_properties =
            mesh->getProperties().createNewPropertyVector<int>(
                field_int_prop_name, MeshLib::MeshItemType::IntegrationPoint);
        field_int_properties->resize(mesh->getNumberOfElements() * 2);
        std::iota(field_int_properties->begin(), field_int_properties->end(),
                  1);

        std::string const point_unsigned_prop_name("PointUnsignedProperty");
        auto point_unsigned_properties =
            mesh->getProperties().createNewPropertyVector<unsigned>(
                point_unsigned_prop_name, MeshLib::MeshItemType::Node,
                mesh->getNumberOfNodes(), 1);
        std::iota(point_unsigned_properties->begin(),
                  point_unsigned_properties->end(),
                  1);

        std::string const cell_unsigned_prop_name("CellUnsignedProperty");
        auto cell_unsigned_properties =
            mesh->getProperties().createNewPropertyVector<unsigned>(
                cell_unsigned_prop_name, MeshLib::MeshItemType::Cell,
                mesh->getNumberOfElements(), 1);
        std::iota(cell_unsigned_properties->begin(),
                  cell_unsigned_properties->end(),
                  1);

        std::string const field_unsigned_prop_name("FieldUnsignedProperty");
        auto field_unsigned_properties =
            mesh->getProperties().createNewPropertyVector<unsigned>(
                field_unsigned_prop_name,
                MeshLib::MeshItemType::IntegrationPoint);
        field_unsigned_properties->resize(mesh->getNumberOfElements() * 2);
        std::iota(field_unsigned_properties->begin(),
                  field_unsigned_properties->end(),
                  1);

        std::string const material_ids_name("MaterialIDs");
        auto material_id_properties =
            mesh->getProperties().createNewPropertyVector<int>(
                material_ids_name, MeshLib::MeshItemType::Cell,
                mesh->getNumberOfElements(), 1);
        std::iota(
            material_id_properties->begin(), material_id_properties->end(), 1);
    }

    ~InSituMesh() override { delete mesh; }
    MeshLib::Mesh* mesh{nullptr};
    const std::size_t subdivisions = 5;
    const double length = 1.0;
    const double dx = length / subdivisions;
};

TEST_F(InSituMesh, Construction)
{
    ASSERT_TRUE(mesh != nullptr);
    ASSERT_EQ((subdivisions + 1) * (subdivisions + 1) * (subdivisions + 1),
              mesh->getNumberOfNodes());
}

// Writes the mesh into a vtk file, reads the file back and converts it into a
// OGS mesh
#ifndef USE_PETSC
TEST_F(InSituMesh, MappedMeshSourceRoundtrip)
#else
TEST_F(InSituMesh, DISABLED_MappedMeshSourceRoundtrip)
#endif
{
    // TODO Add more comparison criteria

    ASSERT_TRUE(mesh != nullptr);
    std::string test_data_file =
        (std::filesystem::temp_directory_path() /= BaseLib::randomString(32))
            .string();

    // -- Test VtkMappedMeshSource, i.e. OGS mesh to VTK mesh
    vtkNew<MeshLib::VtkMappedMeshSource> vtkSource;
    vtkSource->SetMesh(mesh);
    vtkSource->Update();
    vtkUnstructuredGrid* output = vtkSource->GetOutput();

    // Point and cell numbers
    ASSERT_EQ((subdivisions + 1) * (subdivisions + 1) * (subdivisions + 1),
              output->GetNumberOfPoints());
    ASSERT_EQ(subdivisions * subdivisions * subdivisions,
              output->GetNumberOfCells());

    // Point data arrays
    vtkDataArray* pointDoubleArray =
        output->GetPointData()->GetScalars("PointDoubleProperty");
    ASSERT_EQ(pointDoubleArray->GetSize(), mesh->getNumberOfNodes());
    ASSERT_EQ(pointDoubleArray->GetComponent(0, 0), 1.0);
    double* range = pointDoubleArray->GetRange(0);
    ASSERT_EQ(range[0], 1.0);
    ASSERT_EQ(range[1], 1.0 + mesh->getNumberOfNodes() - 1.0);

    vtkDataArray* pointIntArray =
        output->GetPointData()->GetScalars("PointIntProperty");
    ASSERT_EQ(pointIntArray->GetSize(), mesh->getNumberOfNodes());
    ASSERT_EQ(pointIntArray->GetComponent(0, 0), 1.0);
    range = pointIntArray->GetRange(0);
    ASSERT_EQ(range[0], 1.0);
    ASSERT_EQ(range[1], 1 + mesh->getNumberOfNodes() - 1);

    vtkDataArray* pointUnsignedArray =
        output->GetPointData()->GetScalars("PointUnsignedProperty");
    ASSERT_EQ(pointUnsignedArray->GetSize(), mesh->getNumberOfNodes());
    ASSERT_EQ(pointUnsignedArray->GetComponent(0, 0), 1.0);
    range = pointUnsignedArray->GetRange(0);
    ASSERT_EQ(range[0], 1.0);
    ASSERT_EQ(range[1], 1 + mesh->getNumberOfNodes() - 1);

    // Cell data arrays
    vtkDataArray* cellDoubleArray =
        output->GetCellData()->GetScalars("CellDoubleProperty");
    ASSERT_EQ(cellDoubleArray->GetSize(), mesh->getNumberOfElements());
    ASSERT_EQ(cellDoubleArray->GetComponent(0, 0), 1.0);
    range = cellDoubleArray->GetRange(0);
    ASSERT_EQ(range[0], 1.0);
    ASSERT_EQ(range[1], 1.0 + mesh->getNumberOfElements() - 1.0);

    vtkDataArray* cellIntArray =
        output->GetCellData()->GetScalars("CellIntProperty");
    ASSERT_EQ(cellIntArray->GetSize(), mesh->getNumberOfElements());
    ASSERT_EQ(cellIntArray->GetComponent(0, 0), 1.0);
    range = cellIntArray->GetRange(0);
    ASSERT_EQ(range[0], 1.0);
    ASSERT_EQ(range[1], 1 + mesh->getNumberOfElements() - 1);

    vtkDataArray* cellUnsignedArray =
        output->GetCellData()->GetScalars("CellUnsignedProperty");
    ASSERT_EQ(cellUnsignedArray->GetSize(), mesh->getNumberOfElements());
    ASSERT_EQ(cellUnsignedArray->GetComponent(0, 0), 1.0);
    range = cellUnsignedArray->GetRange(0);
    ASSERT_EQ(range[0], 1.0);
    ASSERT_EQ(range[1], 1 + mesh->getNumberOfElements() - 1);

    // Field data arrays
    vtkDataArray* fieldDoubleArray = vtkDataArray::SafeDownCast(
        output->GetFieldData()->GetAbstractArray("FieldDoubleProperty"));
    ASSERT_EQ(fieldDoubleArray->GetSize(), mesh->getNumberOfElements() * 2);
    ASSERT_EQ(fieldDoubleArray->GetComponent(0, 0), 1.0);
    range = fieldDoubleArray->GetRange(0);
    ASSERT_EQ(range[0], 1.0);
    ASSERT_EQ(range[1], mesh->getNumberOfElements() * 2);

    vtkDataArray* fieldIntArray = vtkDataArray::SafeDownCast(
        output->GetFieldData()->GetAbstractArray("FieldIntProperty"));
    ASSERT_EQ(fieldIntArray->GetSize(), mesh->getNumberOfElements() * 2);
    ASSERT_EQ(fieldIntArray->GetComponent(0, 0), 1.0);
    range = fieldIntArray->GetRange(0);
    ASSERT_EQ(range[0], 1.0);
    ASSERT_EQ(range[1], mesh->getNumberOfElements() * 2);

    vtkDataArray* fieldUnsignedArray = vtkDataArray::SafeDownCast(
        output->GetFieldData()->GetAbstractArray("FieldUnsignedProperty"));
    ASSERT_EQ(fieldUnsignedArray->GetSize(), mesh->getNumberOfElements() * 2);
    ASSERT_EQ(fieldUnsignedArray->GetComponent(0, 0), 1.0);
    range = fieldUnsignedArray->GetRange(0);
    ASSERT_EQ(range[0], 1.0);
    ASSERT_EQ(range[1], mesh->getNumberOfElements() * 2);

    // -- Write VTK mesh to file (in all combinations of binary, appended and
    // compressed)
    // TODO: atm vtkXMLWriter::Ascii fails
    for (int dataMode : {vtkXMLWriter::Appended, vtkXMLWriter::Binary})
    {
        for (bool compressed : {true, false})
        {
            if (dataMode == vtkXMLWriter::Ascii && compressed)
            {
                continue;
            }
            MeshLib::IO::VtuInterface vtuInterface(mesh, dataMode, compressed);
            ASSERT_TRUE(vtuInterface.writeToFile(test_data_file));

            // -- Read back VTK mesh
            vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
                vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(test_data_file.c_str());
            reader->Update();
            vtkUnstructuredGrid* vtkMesh = reader->GetOutput();

            // Both VTK meshes should be identical
            ASSERT_EQ(vtkMesh->GetNumberOfPoints(),
                      output->GetNumberOfPoints());
            ASSERT_EQ(vtkMesh->GetNumberOfCells(), output->GetNumberOfCells());
            ASSERT_EQ(vtkMesh->GetFieldData()->GetNumberOfArrays(),
                      output->GetFieldData()->GetNumberOfArrays());
            ASSERT_EQ(vtkMesh->GetPointData()
                          ->GetScalars("PointDoubleProperty")
                          ->GetNumberOfTuples(),
                      pointDoubleArray->GetNumberOfTuples());
            ASSERT_EQ(vtkMesh->GetPointData()
                          ->GetScalars("PointIntProperty")
                          ->GetNumberOfTuples(),
                      pointIntArray->GetNumberOfTuples());
            ASSERT_EQ(vtkMesh->GetPointData()
                          ->GetScalars("PointUnsignedProperty")
                          ->GetNumberOfTuples(),
                      pointUnsignedArray->GetNumberOfTuples());
            ASSERT_EQ(vtkMesh->GetCellData()
                          ->GetScalars("CellDoubleProperty")
                          ->GetNumberOfTuples(),
                      cellDoubleArray->GetNumberOfTuples());
            ASSERT_EQ(vtkMesh->GetCellData()
                          ->GetScalars("CellIntProperty")
                          ->GetNumberOfTuples(),
                      cellIntArray->GetNumberOfTuples());
            ASSERT_EQ(vtkMesh->GetCellData()
                          ->GetScalars("CellUnsignedProperty")
                          ->GetNumberOfTuples(),
                      cellUnsignedArray->GetNumberOfTuples());

            auto get_field_data =
                [&vtkMesh](std::string const array_name) -> vtkDataArray*
            {
                return vtkDataArray::SafeDownCast(
                    vtkMesh->GetFieldData()->GetAbstractArray(
                        array_name.c_str()));
            };

            ASSERT_EQ(
                get_field_data("FieldDoubleProperty")->GetNumberOfTuples(),
                fieldDoubleArray->GetNumberOfTuples());
            ASSERT_EQ(get_field_data("FieldIntProperty")->GetNumberOfTuples(),
                      fieldIntArray->GetNumberOfTuples());
            ASSERT_EQ(
                get_field_data("FieldUnsignedProperty")->GetNumberOfTuples(),
                fieldUnsignedArray->GetNumberOfTuples());

            // Both OGS meshes should be identical
            auto newMesh = std::unique_ptr<MeshLib::Mesh>{
                MeshLib::VtkMeshConverter::convertUnstructuredGrid(vtkMesh)};
            ASSERT_EQ(mesh->getNumberOfNodes(), newMesh->getNumberOfNodes());
            ASSERT_EQ(mesh->getNumberOfElements(),
                      newMesh->getNumberOfElements());

            // Both properties should be identical
            auto meshProperties = mesh->getProperties();
            auto newMeshProperties = newMesh->getProperties();
            ASSERT_EQ(
                newMeshProperties.hasPropertyVector("PointDoubleProperty"),
                meshProperties.hasPropertyVector("PointDoubleProperty"));
            ASSERT_EQ(newMeshProperties.hasPropertyVector("PointIntProperty"),
                      meshProperties.hasPropertyVector("PointIntProperty"));
            ASSERT_EQ(
                newMeshProperties.hasPropertyVector("PointUnsignedProperty"),
                meshProperties.hasPropertyVector("PointUnsignedProperty"));
            ASSERT_EQ(newMeshProperties.hasPropertyVector("CellDoubleProperty"),
                      meshProperties.hasPropertyVector("CellDoubleProperty"));
            ASSERT_EQ(newMeshProperties.hasPropertyVector("CellIntProperty"),
                      meshProperties.hasPropertyVector("CellIntProperty"));
            ASSERT_EQ(
                newMeshProperties.hasPropertyVector("CellUnsignedProperty"),
                meshProperties.hasPropertyVector("CellUnsignedProperty"));

            ASSERT_EQ(
                newMeshProperties.hasPropertyVector("FieldDoubleProperty"),
                meshProperties.hasPropertyVector("FieldDoubleProperty"));
            ASSERT_EQ(newMeshProperties.hasPropertyVector("FieldIntProperty"),
                      meshProperties.hasPropertyVector("FieldIntProperty"));
            ASSERT_EQ(
                newMeshProperties.hasPropertyVector("FieldUnsignedProperty"),
                meshProperties.hasPropertyVector("FieldUnsignedProperty"));

            ASSERT_EQ(newMeshProperties.hasPropertyVector("MaterialIDs"),
                      meshProperties.hasPropertyVector("MaterialIDs"));

            // Check some properties on equality
            auto const* const doubleProps =
                meshProperties.getPropertyVector<double>("PointDoubleProperty");
            auto const* const newDoubleProps =
                newMeshProperties.getPropertyVector<double>(
                    "PointDoubleProperty");
            ASSERT_EQ(newDoubleProps->getNumberOfGlobalComponents(),
                      doubleProps->getNumberOfGlobalComponents());
            ASSERT_EQ(newDoubleProps->getNumberOfTuples(),
                      doubleProps->getNumberOfTuples());
            ASSERT_EQ(newDoubleProps->size(), doubleProps->size());
            for (std::size_t i = 0; i < doubleProps->size(); i++)
            {
                ASSERT_EQ((*newDoubleProps)[i], (*doubleProps)[i]);
            }

            auto unsignedProps = meshProperties.getPropertyVector<unsigned>(
                "CellUnsignedProperty");
            auto newUnsignedIds = newMeshProperties.getPropertyVector<unsigned>(
                "CellUnsignedProperty");

            ASSERT_EQ(newUnsignedIds->getNumberOfGlobalComponents(),
                      unsignedProps->getNumberOfGlobalComponents());
            ASSERT_EQ(newUnsignedIds->getNumberOfTuples(),
                      unsignedProps->getNumberOfTuples());
            ASSERT_EQ(newUnsignedIds->size(), unsignedProps->size());
            for (std::size_t i = 0; i < unsignedProps->size(); i++)
            {
                ASSERT_EQ((*newUnsignedIds)[i], (*unsignedProps)[i]);
            }

            {  // Field data
                auto p =
                    meshProperties.getPropertyVector<int>("FieldIntProperty");
                auto new_p = newMeshProperties.getPropertyVector<int>(
                    "FieldIntProperty");

                ASSERT_EQ(new_p->getNumberOfGlobalComponents(),
                          p->getNumberOfGlobalComponents());
                ASSERT_EQ(new_p->getNumberOfTuples(), p->getNumberOfTuples());
                ASSERT_EQ(new_p->size(), p->size());
                for (std::size_t i = 0; i < unsignedProps->size(); i++)
                {
                    ASSERT_EQ((*newUnsignedIds)[i], (*unsignedProps)[i]);
                }
            }

            auto const* const materialIds =
                meshProperties.getPropertyVector<int>("MaterialIDs");
            auto const* const newMaterialIds =
                newMeshProperties.getPropertyVector<int>("MaterialIDs");
            ASSERT_EQ(newMaterialIds->getNumberOfGlobalComponents(),
                      materialIds->getNumberOfGlobalComponents());
            ASSERT_EQ(newMaterialIds->getNumberOfTuples(),
                      materialIds->getNumberOfTuples());
            ASSERT_EQ(newMaterialIds->size(), materialIds->size());
            for (std::size_t i = 0; i < materialIds->size(); i++)
            {
                ASSERT_EQ((*newMaterialIds)[i], (*materialIds)[i]);
            }
        }
    }
}
