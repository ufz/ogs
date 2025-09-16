/**
 * \file
 * \date 2023-04-25
 * \brief Creates a boundary representation for a layered 3D mesh.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <algorithm>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include "Applications/FileIO/TetGenInterface.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/IO/readStringListFromFile.h"
#include "BaseLib/MPI.h"
#include "GeoLib/IO/AsciiRasterInterface.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshGenerators/LayeredVolume.h"
#include "MeshToolsLib/MeshGenerators/MeshLayerMapper.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Creates a boundary representation for a layered 3D mesh in "
        "*.smesh-format. This boundary representation can be used with the "
        "TetGen Tetrahedral Mesh Generator to create 3D meshes of the given "
        "geometry at arbitrary resolutions and with varying properties. "
        "Details on command line switches and possible parametrisation can be "
        "found in the TetGen User's Manual. Supported raster formats are "
        "ArcGIS ascii rasters (*.asc), Surfer Grids (*.grd) and XYZ raster "
        "files (*.xyz)."
        "Only input meshes consisting of line and triangle elements are "
        "currently supported as mapping of quads might result in invalid mesh "
        "elements.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::SwitchArg use_ascii_arg("", "ascii_output",
                                   "Write VTU output in ASCII format.");
    cmd.add(use_ascii_arg);

    double min_thickness(std::numeric_limits<double>::epsilon());
    TCLAP::ValueArg<double> min_thickness_arg(
        "t", "thickness",
        "The minimum thickness of a layer to be integrated at any given "
        "location (min = 0)",
        false, min_thickness, "MIN_THICKNESS");
    cmd.add(min_thickness_arg);

    TCLAP::ValueArg<std::string> raster_path_arg(
        "r", "raster-list",
        "Input (.vtu). An ascii-file containing a list of raster files, "
        "starting from top "
        "(DEM) to bottom.",
        true, "", "INPUT_FILE_LIST");
    cmd.add(raster_path_arg);

    TCLAP::ValueArg<std::string> mesh_out_arg(
        "o", "output-mesh-file",
        "Output (.smesh). The file name of the resulting 3D mesh", true, "",
        "OUTPUT_FILE");
    cmd.add(mesh_out_arg);

    TCLAP::ValueArg<std::string> mesh_arg(
        "i", "input-mesh-file",
        "Input (.vtu | .msh). The file name of the 2D input mesh", true, "",
        "INPUT_FILE");
    cmd.add(mesh_arg);

    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);

    if (min_thickness_arg.isSet())
    {
        min_thickness = min_thickness_arg.getValue();
        if (min_thickness < 0)
        {
            ERR("Minimum layer thickness must be non-negative value.");
            return EXIT_FAILURE;
        }
    }

    INFO("Reading mesh '{:s}' ... ", mesh_arg.getValue());
    std::unique_ptr<MeshLib::Mesh> const sfc_mesh(MeshLib::IO::readMeshFromFile(
        mesh_arg.getValue(), true /* compute_element_neighbors */));
    if (!sfc_mesh)
    {
        ERR("Error reading mesh '{:s}'.", mesh_arg.getValue());
        return EXIT_FAILURE;
    }
    if (sfc_mesh->getDimension() != 2)
    {
        ERR("Input mesh must be a 2D mesh.");
        return EXIT_FAILURE;
    }
    INFO("done.");

    std::vector<std::string> raster_paths =
        BaseLib::IO::readStringListFromFile(raster_path_arg.getValue());
    if (raster_paths.size() < 2)
    {
        ERR("At least two raster files needed to create 3D mesh.");
        return EXIT_FAILURE;
    }
    std::reverse(raster_paths.begin(), raster_paths.end());

    std::string output_name(mesh_out_arg.getValue());
    if (!BaseLib::hasFileExtension(".smesh", output_name))
    {
        output_name.append(".smesh");
    }

    auto const rasters = FileIO::readRasters(raster_paths);
    LayeredVolume lv;
    if (rasters)
    {
        if (!lv.createLayers(*sfc_mesh, *rasters, min_thickness))
        {
            ERR("Creating the layers was erroneous.");
            return EXIT_FAILURE;
        }
    }
    else
    {
        ERR("The raster files are not accessible.");
        return EXIT_FAILURE;
    }
    std::unique_ptr<MeshLib::Mesh> tg_mesh =
        lv.getMesh("BoundaryRepresentation");
    if (tg_mesh != nullptr)
    {
        std::vector<MeshLib::Node> tg_attr(lv.getAttributePoints());
        FileIO::TetGenInterface tetgen_interface;
        tetgen_interface.writeTetGenSmesh(output_name, *tg_mesh, tg_attr);
        INFO("Smesh was successfully written.");
    }
    else
    {
        ERR("The tetgen-smesh could not be created.");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
