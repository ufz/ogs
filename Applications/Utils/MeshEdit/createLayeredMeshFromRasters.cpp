/*
 * \date 2016-02-11
 * \brief Creates a layered mesh from a 2D mesh and a bunch of raster files.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <algorithm>
#include <fstream>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/FileTools.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "Applications/FileIO/AsciiRasterInterface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshLayerMapper.h"

int readRasterPaths(std::string const& raster_list_file, std::vector<std::string> &raster_path_vec)
{
    std::ifstream in (raster_list_file.c_str());
    if (in.fail())
    {
        ERR ("Could not open file %s.", raster_list_file.c_str());
        return -1;
    }
    std::string line;
    while (getline(in, line))
    {
        if (line.empty())
            continue;
        raster_path_vec.push_back(line);
    }
    if (raster_path_vec.size()<2)
    {
        ERR ("At least two raster files needed to create 3D mesh.");
        return -2;
    }
    std::reverse(raster_path_vec.begin(), raster_path_vec.end());
    return 0;
}

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Creates a layered 3D OGS mesh from an existing 2D OGS mesh and raster "
        "files representing subsurface layers. Supported raster formats are "
        "ArcGIS ascii rasters (*.asc) and Surfer Grids (*.grd)."
        "",
        ' ',
        "1.0");

    TCLAP::ValueArg<std::string> mesh_arg("i", "input-mesh-file",
        "The name of the file containing the 2D input mesh.", true, "", "input file name");
    cmd.add(mesh_arg);

    TCLAP::ValueArg<std::string> mesh_out_arg("o", "output-mesh-file",
        "The name of the file to which the resulting 3D mesh will be written.",
        true, "", "output file name");
    cmd.add(mesh_out_arg);

    TCLAP::ValueArg<std::string> raster_path_arg("r", "raster-list", 
        "An ascii-file containing a list of raster files, starting from top (DEM) to bottom.",
        true, "", "list of raster files");
    cmd.add(raster_path_arg);

    double min_thickness (std::numeric_limits<double>::epsilon());
    TCLAP::ValueArg<double> min_thickness_arg("t", "thickness",
        "The minimum thickness of a layer to be integrated at any given location.",
        false, min_thickness, "minimum layer thickness");
    cmd.add(min_thickness_arg);

    cmd.parse(argc, argv);

    if (min_thickness_arg.isSet())
    {
        min_thickness = min_thickness_arg.getValue();
        if (min_thickness < 0)
        {
            ERR("Minimum layer thickness must be non-negative value.");
            return EXIT_FAILURE;
        }
    }

    INFO("Reading mesh \"%s\" ... ", mesh_arg.getValue().c_str());
    std::unique_ptr<MeshLib::Mesh> const sfc_mesh (MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));
    if (!sfc_mesh) {
        ERR("Error reading mesh \"%s\".", mesh_arg.getValue().c_str());
        return EXIT_FAILURE;
    }
    if (sfc_mesh->getDimension() != 2) {
        ERR("Input mesh needs to be a 2D mesh.");
        return EXIT_FAILURE;
    }
    INFO("done.");

    std::vector<std::string> raster_paths;
    if (readRasterPaths(raster_path_arg.getValue(), raster_paths) != 0)
        return EXIT_FAILURE;

    MeshLib::MeshLayerMapper mapper;
    if (auto rasters = FileIO::readRasters(raster_paths))
    {
        if (!mapper.createLayers(*sfc_mesh, *rasters, min_thickness))
            return EXIT_FAILURE;
    }
    else
        return EXIT_FAILURE;

    std::string output_name (mesh_out_arg.getValue());
    if (!BaseLib::hasFileExtension("vtu", output_name))
        output_name.append(".vtu");
    INFO("Writing mesh \"%s\" ... ", output_name.c_str());
    MeshLib::IO::writeMeshToFile(*(mapper.getMesh("SubsurfaceMesh").release()), output_name);
    INFO("done.");

    return EXIT_SUCCESS;
}
