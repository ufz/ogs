/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <memory>

#include "gtest/gtest.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkUnsignedIntArray.h>
#include <vtkUnstructuredGrid.h>

#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/VtkMeshConverter.h"

/// Constructs this VTK mesh: https://gist.github.com/bilke/8bc7e5e7084ad4d806e1
class TestVtkMeshConverter : public ::testing::Test
{
public:
    TestVtkMeshConverter()
        : vtu(nullptr)
    {
        // Points and cells
        auto points = vtkSmartPointer<vtkPoints>::New();
        points->InsertNextPoint(703.875, 758.75, 0.9971);
        points->InsertNextPoint(767.625, 679.5, 6.2356);
        points->InsertNextPoint(835.25, 602.50, 11.2227);
        points->InsertNextPoint(608.75, 629.375, 10.3723);
        points->InsertNextPoint(672, 549.375, 17.4011);
        points->InsertNextPoint(739.75, 472.5, 20.9258);
        points->InsertNextPoint(703.875, 758.75, 101.07145);
        points->InsertNextPoint(767.625, 679.5, 106.2616);
        points->InsertNextPoint(835.25, 602.5, 111.2233);
        points->InsertNextPoint(608.75, 629.375, 110.4138);
        points->InsertNextPoint(672, 549.375, 117.4165);
        points->InsertNextPoint(739.75, 472.5, 120.92995);

        auto cells = vtkSmartPointer<vtkCellArray>::New();
        auto cell1 = vtkSmartPointer<vtkIdList>::New();
        cell1->InsertNextId(4);
        cell1->InsertNextId(1);
        cell1->InsertNextId(0);
        cell1->InsertNextId(3);
        cell1->InsertNextId(10);
        cell1->InsertNextId(7);
        cell1->InsertNextId(6);
        cell1->InsertNextId(9);
        auto cell2 = vtkSmartPointer<vtkIdList>::New();
        cell2->InsertNextId(5);
        cell2->InsertNextId(2);
        cell2->InsertNextId(1);
        cell2->InsertNextId(4);
        cell2->InsertNextId(11);
        cell2->InsertNextId(8);
        cell2->InsertNextId(7);
        cell2->InsertNextId(10);
        cells->InsertNextCell(cell1);
        cells->InsertNextCell(cell2);

        vtu = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtu->SetPoints(points);
        vtu->SetCells(12, cells);

        // Data arrays
        auto pointIds = vtkSmartPointer<vtkIntArray>::New();
        pointIds->SetNumberOfComponents(1);
        pointIds->SetName("PointIDs");
        for (int i = 0; i < 12; ++i)
            pointIds->InsertNextTuple1(i);
        vtu->GetPointData()->AddArray(pointIds);

        auto materialIds = vtkSmartPointer<vtkUnsignedIntArray>::New();
        materialIds->SetNumberOfComponents(1);
        materialIds->SetName("MaterialIDs");
        materialIds->InsertNextTuple1(0);
        materialIds->InsertNextTuple1(1);
        vtu->GetCellData()->AddArray(materialIds);
    }

    static std::size_t const mesh_size = 5;
    vtkSmartPointer<vtkUnstructuredGrid> vtu;
};

TEST_F(TestVtkMeshConverter, Conversion)
{
    auto mesh = std::unique_ptr<MeshLib::Mesh>(
        MeshLib::VtkMeshConverter::convertUnstructuredGrid(vtu));
    ASSERT_EQ(mesh->getNumberOfNodes(), vtu->GetNumberOfPoints());
    ASSERT_EQ(mesh->getNumberOfElements(), vtu->GetNumberOfCells());
    ASSERT_EQ(mesh->getElement(0)->getCellType(), MeshLib::CellType::HEX8);

    auto meshProperties = mesh->getProperties();

    // MaterialIDs are converted to an int property
    auto const* const materialIds =
        meshProperties.getPropertyVector<int>("MaterialIDs");
    ASSERT_NE(nullptr, materialIds);
    auto vtkMaterialIds = vtu->GetCellData()->GetArray("MaterialIDs");
    ASSERT_EQ(materialIds->size(), vtkMaterialIds->GetNumberOfTuples());
    for(std::size_t i = 0; i < materialIds->size(); i++)
        ASSERT_EQ((*materialIds)[i], vtkMaterialIds->GetTuple1(i));
}
