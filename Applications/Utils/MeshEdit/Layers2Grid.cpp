/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "BaseLib/IO/readStringListFromFile.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
#include "GeoLib/AABB.h"
#include "InfoLib/GitInfo.h"
#include "MathLib/Point3d.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshToolsLib/MeshGenerators/VoxelGridFromLayeredMeshes.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Reads a list of 2D unstructured mesh layers and samples them onto a "
        "structured grid of the same extent. Note, that a large cube size may "
        "result in an undersampling of the original structure.\nCube sizes are "
        "defines by x/y/z-parameters. For equilateral cubes, only the "
        "x-parameter needs to be set.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::SwitchArg dilate_arg(
        "d", "dilate",
        "assign mat IDs based on single nodes instead of a majority of nodes, "
        "which can result in a slightly increased voxel grid extent",
        false);
    cmd.add(dilate_arg);

    TCLAP::ValueArg<double> z_arg(
        "z", "cellsize-z",
        "edge length of cubes in z-direction (depth), "
        "(min = 0)",
        false, 1000, "CELLSIZE-Z");
    cmd.add(z_arg);

    TCLAP::ValueArg<double> y_arg(
        "y", "cellsize-y",
        "edge length of cubes in y-direction (latitude), "
        "(min = 0)",
        false, 1000, "CELLSIZE-Y");
    cmd.add(y_arg);

    TCLAP::ValueArg<double> x_arg(
        "x", "cellsize-x",
        "edge length of cubes in x-direction (longitude) or all directions, if "
        "y and z are not set, (min = 0)",
        true, 1000, "CELLSIZE-X");
    cmd.add(x_arg);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "Output (.vtu). Name of output mesh file", true, "",
        "OUTPUT_FILE");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> input_arg(
        "i", "input",
        "Input (.vtu). Name of the input file list containing the paths the "
        "all input layers "
        "in correct order from top to bottom",
        true, "", "INPUT_FILE_LIST");
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
        return EXIT_FAILURE;
    }

    double const x_size = x_arg.getValue();
    double const y_size = (y_arg.isSet()) ? y_arg.getValue() : x_arg.getValue();
    double const z_size = (z_arg.isSet()) ? z_arg.getValue() : x_arg.getValue();
    std::array<double, 3> const cellsize = {x_size, y_size, z_size};

    std::string const input_name = input_arg.getValue();
    std::string const output_name = output_arg.getValue();
    auto const layer_names = BaseLib::IO::readStringListFromFile(input_name);
    if (layer_names.size() < 2)
    {
        ERR("At least two layers are required to create a 3D Mesh");
        return EXIT_FAILURE;
    }

    std::vector<std::unique_ptr<MeshLib::Mesh>> layers;
    layers.reserve(layer_names.size());
    constexpr double minval = std::numeric_limits<double>::max();
    constexpr double maxval = std::numeric_limits<double>::lowest();
    std::pair<MathLib::Point3d, MathLib::Point3d> extent(
        MathLib::Point3d{{minval, minval, minval}},
        MathLib::Point3d{{maxval, maxval, maxval}});

    for (auto const& layer : layer_names)
    {
        auto mesh(MeshLib::IO::readMeshFromFile(layer));
        if (mesh == nullptr)
        {
            ERR("Input layer '{:s}' not found. Aborting...", layer);
            return EXIT_FAILURE;
        }
        layers.emplace_back(mesh);
    }
    std::vector<MeshLib::Mesh const*> layers_ptr;
    std::transform(std::begin(layers), std::end(layers),
                   std::back_inserter(layers_ptr),
                   [](auto const& layer) { return layer.get(); });
    bool const dilate = dilate_arg.getValue();
    auto mesh = MeshToolsLib::MeshGenerators::VoxelFromLayeredMeshes::
        createVoxelFromLayeredMesh(extent, layers_ptr, cellsize, dilate);
    if (mesh == nullptr)
    {
        ERR("The VoxelGrid could not be created.");
        return EXIT_FAILURE;
    }
    MeshLib::IO::VtuInterface vtu(mesh.get());
    vtu.writeToFile(output_name);
    return EXIT_SUCCESS;
}
