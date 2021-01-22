/*
 * \date 2016-02-11
 * \brief Creates a layered mesh from a 2D mesh and a bunch of raster files.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <algorithm>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include <tclap/CmdLine.h>

#include "InfoLib/GitInfo.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/IO/readStringListFromFile.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "Applications/FileIO/AsciiRasterInterface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshGenerators/MeshLayerMapper.h"

int main (int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Creates a layered 3D OGS mesh from an existing 2D OGS mesh and a list "
        "of raster files representing subsurface layers. Supported raster "
        "formats are ArcGIS ascii rasters (*.asc) and Surfer Grids (*.grd). "
        "Only input meshes consisting of line and triangle elements are "
        "currently supported as mapping of quads might result in invalid mesh "
        "elements.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::SwitchArg use_ascii_arg("", "ascii_output",
                                   "Write VTU output in ASCII format.");
    cmd.add(use_ascii_arg);

    double min_thickness(std::numeric_limits<double>::epsilon());
    TCLAP::ValueArg<double> min_thickness_arg(
        "t", "thickness",
        "The minimum thickness of a layer to be integrated at any given "
        "location.",
        false, min_thickness, "floating point number");
    cmd.add(min_thickness_arg);

    TCLAP::ValueArg<std::string> raster_path_arg(
        "r", "raster-list",
        "An ascii-file containing a list of raster files, starting from top "
        "(DEM) to bottom.",
        true, "", "file name");
    cmd.add(raster_path_arg);

        TCLAP::ValueArg<std::string> mesh_out_arg(
        "o", "output-mesh-file",
        "The file name of the resulting 3D mesh.",
        true, "", "file name");
    cmd.add(mesh_out_arg);

    TCLAP::ValueArg<std::string> mesh_arg(
        "i", "input-mesh-file",
        "The file name of the 2D input mesh.", true, "",
        "file name");
    cmd.add(mesh_arg);

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

    INFO("Reading mesh '{:s}' ... ", mesh_arg.getValue());
    std::unique_ptr<MeshLib::Mesh> const sfc_mesh(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));
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
    if (raster_paths.size()<2)
    {
        ERR ("At least two raster files needed to create 3D mesh.");
        return EXIT_FAILURE;
    }
    std::reverse(raster_paths.begin(), raster_paths.end());

    MeshLib::MeshLayerMapper mapper;
    if (auto rasters = FileIO::readRasters(raster_paths))
    {
        if (!mapper.createLayers(*sfc_mesh, *rasters, min_thickness))
        {
            return EXIT_FAILURE;
        }
    }
    else
    {
        return EXIT_FAILURE;
    }

    std::string output_name(mesh_out_arg.getValue());
    if (!BaseLib::hasFileExtension(".vtu", output_name))
    {
        output_name.append(".vtu");
    }

    INFO("Writing mesh '{:s}' ... ", output_name);
    auto const result_mesh = mapper.getMesh("SubsurfaceMesh");
    if (result_mesh == nullptr)
    {
        ERR("Mapper returned empty result for 'SubsurfaceMesh'.");
        return EXIT_FAILURE;
    }

    auto const data_mode =
        use_ascii_arg.getValue() ? vtkXMLWriter::Ascii : vtkXMLWriter::Binary;

    MeshLib::IO::writeVtu(*result_mesh, output_name, data_mode);
    INFO("done.");

    return EXIT_SUCCESS;
}
