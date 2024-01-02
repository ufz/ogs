/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#include "VoxelGridFromMesh.h"

#include <vtkAbstractArray.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkLongLongArray.h>
#include <vtkXMLUnstructuredGridReader.h>

#include <cassert>
#include <cmath>

#include "InfoLib/GitInfo.h"
#include "MathLib/Point3d.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/MeshSearch/MeshElementGrid.h"
#include "MeshToolsLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"

namespace MeshToolsLib::MeshGenerator::VoxelGridFromMesh
{
std::array<std::size_t, 3> getNumberOfVoxelPerDimension(
    std::array<double, 3> const& ranges, std::array<double, 3> const& cellsize)
{
    if (cellsize[0] <= 0 || cellsize[1] <= 0 || cellsize[2] <= 0)
    {
        OGS_FATAL("A cellsize ({},{},{}) is not allowed to be <= 0",
                  cellsize[0], cellsize[1], cellsize[2]);
    }
    std::array<double, 3> numberOfVoxel = {ranges[0] / cellsize[0],
                                           ranges[1] / cellsize[1],
                                           ranges[2] / cellsize[2]};

    if (ranges[0] < 0 || ranges[1] < 0 || ranges[2] < 0)
    {
        OGS_FATAL(
            "The difference of max-min ({},{},{}) is not allowed to be < 0",
            ranges[0], ranges[1], ranges[2]);
    }
    std::replace(numberOfVoxel.begin(), numberOfVoxel.end(), 0, 1);

    return {static_cast<std::size_t>(std::lround(numberOfVoxel[0])),
            static_cast<std::size_t>(std::lround(numberOfVoxel[1])),
            static_cast<std::size_t>(std::lround(numberOfVoxel[2]))};
}

std::vector<int> assignCellIds(vtkSmartPointer<vtkUnstructuredGrid> const& mesh,
                               MathLib::Point3d const& min,
                               std::array<std::size_t, 3> const& dims,
                               std::array<double, 3> const& cellsize)
{
    vtkSmartPointer<vtkCellLocator> locator =
        vtkSmartPointer<vtkCellLocator>::New();
    locator->SetDataSet(mesh);
    locator->Update();

    std::vector<int> cell_ids;
    cell_ids.reserve(dims[0] * dims[1] * dims[2]);
    std::array<double, 3> const grid_max = {min[0] + dims[0] * cellsize[0],
                                            min[1] + dims[1] * cellsize[1],
                                            min[2] + dims[2] * cellsize[2]};

    double const start[3] = {min[0] + cellsize[0] / 2.,
                             min[1] + cellsize[1] / 2.,
                             min[2] + cellsize[2] / 2.};
    double pnt[3];
    for (pnt[2] = start[2]; pnt[2] < grid_max[2]; pnt[2] += cellsize[2])
    {
        for (pnt[1] = start[1]; pnt[1] < grid_max[1]; pnt[1] += cellsize[1])
        {
            for (pnt[0] = start[0]; pnt[0] < grid_max[0]; pnt[0] += cellsize[0])
            {
                cell_ids.push_back(static_cast<int>(locator->FindCell(pnt)));
            }
        }
    }
    return cell_ids;
}

bool removeUnusedGridCells(vtkSmartPointer<vtkUnstructuredGrid> const& mesh,
                           std::unique_ptr<MeshLib::Mesh>& grid)
{
    MeshLib::ElementSearch search(*grid);
    std::size_t const n_elems_marked = search.searchByPropertyValueRange<int>(
        cell_id_name, 0, static_cast<int>(mesh->GetNumberOfCells()), true);

    if (n_elems_marked == grid->getNumberOfElements())
    {
        ERR("No valid elements found. Aborting...");
        return false;
    }

    if (n_elems_marked)
    {
        grid.reset(MeshToolsLib::removeElements(
            *grid, search.getSearchedElementIDs(), "trimmed_grid"));
    }
    return true;
}

// PropertyVector is required to contain parameter cell_id_name
template <typename T, typename VTK_TYPE>
void mapArray(MeshLib::Mesh& grid, VTK_TYPE vtk_arr, std::string const& name)
{
    auto const* cell_ids = grid.getProperties().getPropertyVector<int>(
        cell_id_name, MeshLib::MeshItemType::Cell, 1);
    if (cell_ids == nullptr)
    {
        // Error message
        return;
    }
    auto& arr = *grid.getProperties().createNewPropertyVector<T>(
        name, MeshLib::MeshItemType::Cell, vtk_arr->GetNumberOfComponents());
    std::size_t const n_elems = cell_ids->size();
    arr.resize(n_elems);
    for (std::size_t j = 0; j < n_elems; ++j)
        arr[j] = vtk_arr->GetValue((*cell_ids)[j]);
}

template <typename T>
// check whether dynamic_cast of cell data is possible the type of cell data to
// map them on voxelgrid
bool checkDyncast(MeshLib::Mesh& mesh,
                  vtkSmartPointer<vtkCellData> const cell_data,
                  char const* const name)
{
    using DataArrayType = vtkAOSDataArrayTemplate<T>;
    vtkSmartPointer<DataArrayType> const arr =
        dynamic_cast<DataArrayType*>(cell_data->GetArray(name));
    if (!arr)
    {
        return false;
    }
    mapArray<T, vtkSmartPointer<DataArrayType>>(mesh, arr, name);
    return true;
}

void mapMeshArraysOntoGrid(vtkSmartPointer<vtkUnstructuredGrid> const& mesh,
                           std::unique_ptr<MeshLib::Mesh> const& grid)
{
    assert(mesh != nullptr);
    assert(grid != nullptr);
    vtkSmartPointer<vtkCellData> const cell_data = mesh->GetCellData();
    for (int i = 0; i < cell_data->GetNumberOfArrays(); ++i)
    {
        auto const name = cell_data->GetArrayName(i);

        if (!(checkDyncast<double>(*grid, cell_data, name) ||
              checkDyncast<long>(*grid, cell_data, name) ||
              checkDyncast<long long>(*grid, cell_data, name) ||
              checkDyncast<int>(*grid, cell_data, name)))
        {
            WARN("Ignoring array '{:s}', array type {:s} not implemented...",
                 name,
                 cell_data->GetArray(name)->GetDataTypeAsString());
        }
    }
}
}  // namespace MeshToolsLib::MeshGenerator::VoxelGridFromMesh
