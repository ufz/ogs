/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "VoxelGridFromMesh.h"

namespace MeshLib::MeshGenerator::VoxelGridFromMesh
{

std::array<std::size_t, 3> getDimensions(MathLib::Point3d const& min,
                                         MathLib::Point3d const& max,
                                         std::array<double, 3> const& cellsize)
{
    return {
        static_cast<std::size_t>(std::ceil((max[0] - min[0]) / cellsize[0])),
        static_cast<std::size_t>(std::ceil((max[1] - min[1]) / cellsize[1])),
        static_cast<std::size_t>(std::ceil((max[2] - min[2]) / cellsize[2]))};
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
    for (double k = min[2] + (cellsize[2] / 2.0); k < grid_max[2];
         k += cellsize[2])
    {
        for (double j = min[1] + (cellsize[1] / 2.0); j < grid_max[1];
             j += cellsize[1])
        {
            for (double i = min[0] + (cellsize[0] / 2.0); i < grid_max[0];
                 i += cellsize[0])
            {
                double pnt[3] = {i, j, k};
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
        grid.reset(MeshLib::removeElements(
            *grid, search.getSearchedElementIDs(), "trimmed_grid"));
    }
    return true;
}

template <typename T, typename VTK_TYPE>
void mapArray(MeshLib::Mesh& grid, VTK_TYPE vtk_arr, std::string const& name)
{
    MeshLib::PropertyVector<int> const& cell_ids =
        *grid.getProperties().getPropertyVector<int>(
            cell_id_name, MeshLib::MeshItemType::Cell, 1);
    std::vector<T>& arr = *grid.getProperties().createNewPropertyVector<T>(
        name, MeshLib::MeshItemType::Cell, vtk_arr->GetNumberOfComponents());
    std::size_t const n_elems = cell_ids.size();
    arr.resize(n_elems);
    for (std::size_t j = 0; j < n_elems; ++j)
        arr[j] = vtk_arr->GetValue(cell_ids[j]);
}

void mapMeshArraysOntoGrid(vtkSmartPointer<vtkUnstructuredGrid> const& mesh,
                           std::unique_ptr<MeshLib::Mesh>& grid)
{
    vtkSmartPointer<vtkCellData> const cell_data = mesh->GetCellData();
    for (int i = 0; i < cell_data->GetNumberOfArrays(); ++i)
    {
        std::string const& name = cell_data->GetArrayName(i);
        vtkSmartPointer<vtkDoubleArray> const dbl_arr =
            dynamic_cast<vtkDoubleArray*>(cell_data->GetArray(name.c_str()));
        if (dbl_arr)
        {
            mapArray<double, vtkSmartPointer<vtkDoubleArray>>(*grid, dbl_arr,
                                                              name);
            continue;
        }
        vtkSmartPointer<vtkIntArray> const int_arr =
            dynamic_cast<vtkIntArray*>(cell_data->GetArray(name.c_str()));
        if (int_arr)
        {
            mapArray<int, vtkSmartPointer<vtkIntArray>>(*grid, int_arr, name);
            continue;
        }
        WARN("Ignoring array '{:s}', array type not implemented...", name);
    }
}
}  // namespace MeshLib::MeshGenerator::VoxelGridFromMesh