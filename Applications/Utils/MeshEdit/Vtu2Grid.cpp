/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>
#include <vtkXMLUnstructuredGridReader.h>

#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshToolsLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"
#include "MeshToolsLib/MeshGenerators/VoxelGridFromMesh.h"

#ifdef USE_PETSC
#include <mpi.h>
#endif
int main(int argc, char* argv[])
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
            "Copyright (c) 2012-2024, OpenGeoSys Community "
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

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    if ((y_arg.isSet() && !z_arg.isSet()) ||
        ((!y_arg.isSet() && z_arg.isSet())))
    {
        ERR("For equilateral cubes, only x needs to be set. For unequal "
            "cuboids, all three edge lengths (x/y/z) need to be specified.");
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return -1;
    }
    using namespace MeshToolsLib::MeshGenerator;

    double const x_size = x_arg.getValue();
    double const y_size = (y_arg.isSet()) ? y_arg.getValue() : x_arg.getValue();
    double const z_size = (z_arg.isSet()) ? z_arg.getValue() : x_arg.getValue();

    if (x_size <= 0 || y_size <= 0 || z_size <= 0)
    {
        ERR("A cellsize ({},{},{}) is not allowed to be <= 0", x_size, y_size,
            z_size);
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return -1;
    }

    std::array<double, 3> const cellsize = {x_size, y_size, z_size};

    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader =
        vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(input_arg.getValue().c_str());
    reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> mesh = reader->GetOutput();

    double* const bounds = mesh->GetBounds();
    MathLib::Point3d const min(
        std::array<double, 3>{bounds[0], bounds[2], bounds[4]});
    MathLib::Point3d const max(
        std::array<double, 3>{bounds[1], bounds[3], bounds[5]});
    std::array<double, 3> ranges = {max[0] - min[0], max[1] - min[1],
                                    max[2] - min[2]};
    if (ranges[0] < 0 || ranges[1] < 0 || ranges[2] < 0)
    {
        ERR("The range ({},{},{}) is not allowed to be < 0", ranges[0],
            ranges[1], ranges[2]);
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return -1;
    }
    std::array<std::size_t, 3> const dims =
        VoxelGridFromMesh::getNumberOfVoxelPerDimension(ranges, cellsize);
    std::unique_ptr<MeshLib::Mesh> grid(
        MeshToolsLib::MeshGenerator::generateRegularHexMesh(
            dims[0], dims[1], dims[2], cellsize[0], cellsize[1], cellsize[2],
            min, "grid"));

    std::vector<int> const tmp_ids =
        VoxelGridFromMesh::assignCellIds(mesh, min, dims, cellsize);
    std::vector<int>& cell_ids =
        *grid->getProperties().createNewPropertyVector<int>(
            VoxelGridFromMesh::cell_id_name, MeshLib::MeshItemType::Cell, 1);
    std::copy(tmp_ids.cbegin(), tmp_ids.cend(), std::back_inserter(cell_ids));

    if (!VoxelGridFromMesh::removeUnusedGridCells(mesh, grid))
    {
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

    VoxelGridFromMesh::mapMeshArraysOntoGrid(mesh, grid);

    if (MeshLib::IO::writeMeshToFile(*grid, output_arg.getValue()) != 0)
    {
#ifdef USE_PETSC
        MPI_Finalize();
#endif
        return EXIT_FAILURE;
    }

#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
