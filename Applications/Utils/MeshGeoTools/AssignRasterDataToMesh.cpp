// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <tclap/CmdLine.h>

#include <memory>
#include <string>

#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
#include "GeoLib/IO/AsciiRasterInterface.h"
#include "GeoLib/Raster.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshEditing/RasterDataToMesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Assigns the pixel information of a raster file to a scalar array of a "
        "specified 2D mesh. Data will be assigned to a node array by default. "
        "Adding information to cell arrays is also possible, pixel values at "
        "the centre of the cell will be used in this case. Note that large "
        "differences in resolution between cell size of the mesh and pixel "
        "size of the raster can give unexpected results. A no-data value will "
        "be added in case of missing or transparent values.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2026, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<double> nodata_arg(
        "e", "nodata",
        "The no data value used for missing values "
        "(min = 0)",
        false, 0, "NO_DATA");
    cmd.add(nodata_arg);

    TCLAP::SwitchArg set_cells_arg("c", "cell-array",
                                   "Assigns raster data to cell array");
    cmd.add(set_cells_arg);

    TCLAP::SwitchArg set_nodes_arg(
        "n", "node-array", "Assigns raster data to node array (default)");
    cmd.add(set_nodes_arg);

    TCLAP::ValueArg<std::string> array_name_arg(
        "s", "scalar-name", "The name of the newly created scalar array.", true,
        "", "SCALAR_NAME");
    cmd.add(array_name_arg);
    TCLAP::ValueArg<std::string> raster_arg(
        "r", "raster", "Input (.asc). Name of the input raster file", true, "",
        "INPUT_FILE");
    cmd.add(raster_arg);

    TCLAP::ValueArg<std::string> output_arg(
        "o", "output", "Output (.vtu). Name of the output mesh file", true, "",
        "OUTPUT_FILE");
    cmd.add(output_arg);
    TCLAP::ValueArg<std::string> input_arg(
        "i", "input", "Input (.vtu). Name of the input mesh file", true, "",
        "INPUT_FILE");
    cmd.add(input_arg);
    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    bool const create_cell_array(set_cells_arg.isSet());
    bool const create_node_array =
        (create_cell_array) ? set_nodes_arg.isSet() : true;

    std::string const& mesh_name = input_arg.getValue();
    std::string const& output_name = output_arg.getValue();
    std::string const& raster_name = raster_arg.getValue();

    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshLib::IO::readMeshFromFile(mesh_name));
    if (mesh->getDimension() > 2)
    {
        ERR("Method can not be applied to 3D meshes.");
        return EXIT_FAILURE;
    }

    std::unique_ptr<GeoLib::Raster> const raster(
        FileIO::AsciiRasterInterface::getRasterFromASCFile(raster_name));

    if (create_node_array)
    {
        bool const assigned = MeshToolsLib::RasterDataToMesh::projectToNodes(
            *mesh, *raster, nodata_arg.getValue(), array_name_arg.getValue());

        if (!assigned)
        {
            ERR("Error assigning raster data to scalar node array");
            return EXIT_FAILURE;
        }
        INFO("Created node array {:s}", array_name_arg.getValue());
    }

    if (create_cell_array)
    {
        bool const assigned = MeshToolsLib::RasterDataToMesh::projectToElements(
            *mesh, *raster, nodata_arg.getValue(), array_name_arg.getValue());

        if (!assigned)
        {
            ERR("Error assigning raster data to scalar cell array");
            return EXIT_FAILURE;
        }
        INFO("Created cell array {:s}", array_name_arg.getValue());
    }

    MeshLib::IO::VtuInterface vtu(mesh.get());
    vtu.writeToFile(output_name);
    return EXIT_SUCCESS;
}
