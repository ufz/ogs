/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include "InfoLib/GitInfo.h"

#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshLib/MeshGenerators/MeshGenerator.h"
#include "MeshLib/MeshSearch/ElementSearch.h"

#include <vtkAbstractArray.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>

std::string const cell_id_name = "CellIds";

std::array<std::size_t, 3> getDimensions(MathLib::Point3d const& min,
                                         MathLib::Point3d const& max,
                                         std::array<double, 3> const& cellsize)
{
    return
    {
        static_cast<std::size_t>(std::ceil((max[0] - min[0]) / cellsize[0])),
        static_cast<std::size_t>(std::ceil((max[1] - min[1]) / cellsize[1])),
        static_cast<std::size_t>(std::ceil((max[2] - min[2]) / cellsize[2]))
    };
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

    std::vector<int>  cell_ids;
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
        ERR ("No valid elements found. Aborting...");
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
            mapArray<double, vtkSmartPointer<vtkDoubleArray>>(*grid, dbl_arr, name);
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

int main (int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Reads a 3D unstructured mesh and samples it onto a structured grid of "
        "the same extent. Cell properties are mapped onto the grid (sampled at "
        "the centre-points of each cube), node properties are ignored. Note, "
        "that a large cube size may result in an undersampling of the original "
        "mesh structure.\nCube sizes are defines by x/y/z-parameters. For "
        "equilateral cubes, only the x-parameter needs to be set.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<double> z_arg("z", "cellsize-z",
                                  "edge length of cubes in z-direction (depth)",
                                  false, 1000, "floating point number");
    cmd.add(z_arg);

    TCLAP::ValueArg<double> y_arg(
        "y", "cellsize-y", "edge length of cubes in y-direction (latitude)",
        false, 1000, "floating point number");
    cmd.add(y_arg);

    TCLAP::ValueArg<double> x_arg(
        "x", "cellsize-x",
        "edge length of cubes in x-direction (longitude) or all directions, if "
        "y and z are not set",
        true, 1000, "floating point number");
    cmd.add(x_arg);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "the output grid (*.vtu)", true, "", "output.vtu");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> input_arg("i", "input",
                                           "the 3D input mesh (*.vtu, *.msh)",
                                           true, "", "input.vtu");
    cmd.add(input_arg);
    cmd.parse(argc, argv);

    if ((y_arg.isSet() && !z_arg.isSet()) ||
        ((!y_arg.isSet() && z_arg.isSet())))
    {
        ERR("For equilateral cubes, only x needs to be set. For unequal "
            "cuboids, all three edge lengths (x/y/z) need to be specified.");
        return -1;
    }

    double const x_size = x_arg.getValue();
    double const y_size = (y_arg.isSet()) ? y_arg.getValue() : x_arg.getValue();
    double const z_size = (z_arg.isSet()) ? z_arg.getValue() : x_arg.getValue();
    std::array<double, 3> const cellsize = { x_size, y_size, z_size };

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(input_arg.getValue().c_str());
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();

    double* const bounds = mesh->GetBounds();
    MathLib::Point3d const min(std::array<double, 3>{bounds[0], bounds[2], bounds[4]});
    MathLib::Point3d const max(std::array<double, 3>{bounds[1], bounds[3], bounds[5]});
    std::array<std::size_t, 3> const dims = getDimensions(min, max, cellsize);
    std::unique_ptr<MeshLib::Mesh> grid(
        MeshLib::MeshGenerator::generateRegularHexMesh(
            dims[0], dims[1], dims[2], cellsize[0], cellsize[1], cellsize[2],
            min, "grid"));


    std::vector<int> const tmp_ids = assignCellIds(mesh, min, dims, cellsize);
    std::vector<int>& cell_ids =
        *grid->getProperties().createNewPropertyVector<int>(
            cell_id_name, MeshLib::MeshItemType::Cell, 1);
    std::copy(tmp_ids.cbegin(), tmp_ids.cend(), std::back_inserter(cell_ids));

    if (!removeUnusedGridCells(mesh, grid))
        return EXIT_FAILURE;

    mapMeshArraysOntoGrid(mesh, grid);

    if (MeshLib::IO::writeMeshToFile(*grid, output_arg.getValue()) != 0)
        return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
