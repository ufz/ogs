// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <tclap/CmdLine.h>

#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshToolsLib/MeshEditing/RemoveMeshComponents.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"
#include "MeshToolsLib/MeshGenerators/VoxelGridFromMesh.h"

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
            "Copyright (c) 2012-2026, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<double> z_arg("z", "cellsize-z",
                                  "edge length of cubes in z-direction "
                                  "(depth), (min = 0)",
                                  false, 1000, "CELLSIZE_Z");
    cmd.add(z_arg);

    TCLAP::ValueArg<double> y_arg("y", "cellsize-y",
                                  "edge length of cubes in y-direction "
                                  "(latitude), (min = 0)",
                                  false, 1000, "CELLSIZE_Y");
    cmd.add(y_arg);

    TCLAP::ValueArg<double> x_arg(
        "x", "cellsize-x",
        "edge length of cubes in x-direction (longitude) or all directions, if "
        "y and z are not set, (min = 0)",
        true, 1000, "CELLSIZE_X");
    cmd.add(x_arg);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "Output (.vtu). The output grid file", true, "",
        "OUTPUT_FILE");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> input_arg(
        "i", "input", "Input (.vtu | .msh). The 3D input mesh file", true, "",
        "INPUT_FILE");
    cmd.add(input_arg);
    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    if ((y_arg.isSet() && !z_arg.isSet()) ||
        ((!y_arg.isSet() && z_arg.isSet())))
    {
        ERR("For equilateral cubes, only x needs to be set. For unequal "
            "cuboids, all three edge lengths (x/y/z) need to be specified.");
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
        return -1;
    }

    std::array<double, 3> const cellsize = {x_size, y_size, z_size};

    vtkSmartPointer<vtkUnstructuredGrid> mesh =
        MeshLib::IO::VtuInterface::readVtuFileToVtkUnstructuredGrid(
            input_arg.getValue());
    if (mesh == nullptr)
    {
        return EXIT_FAILURE;
    }

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
    auto* const cell_ids = grid->getProperties().createNewPropertyVector<int>(
        VoxelGridFromMesh::cell_id_name, MeshLib::MeshItemType::Cell, 1);
    assert(cell_ids);
    cell_ids->assign(tmp_ids);

    if (!VoxelGridFromMesh::removeUnusedGridCells(mesh, grid))
    {
        return EXIT_FAILURE;
    }

    VoxelGridFromMesh::mapMeshArraysOntoGrid(mesh, grid);

    if (MeshLib::IO::writeMeshToFile(*grid, output_arg.getValue()) != 0)
    {
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
