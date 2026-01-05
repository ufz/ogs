// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <tclap/CmdLine.h>

#include <algorithm>
#include <iterator>
#include <memory>
#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/join.hpp>
#include <range/v3/view/repeat_n.hpp>
#include <string>
#include <vector>

#include "BaseLib/FileTools.h"
#include "BaseLib/IO/readStringListFromFile.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
#include "GeoLib/IO/AsciiRasterInterface.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshGenerators/MeshLayerMapper.h"

void assignSurfaceMaterialIDsToSubsurfaceLayers(MeshLib::Mesh* sfc_mesh,
                                                MeshLib::Mesh* subsurface_mesh)
{
    auto const* surface_material_ids = materialIDs(*sfc_mesh);
    auto* subsurface_material_ids = materialIDs(*subsurface_mesh);
    // reset the material ids
    if (!surface_material_ids)
    {
        OGS_FATAL("Surface mesh does not contain material IDs");
    }
    if (surface_material_ids->empty())
    {
        OGS_FATAL("Surface mesh material IDs doesn't contain any values");
    }
    if (!subsurface_material_ids)
    {
        OGS_FATAL("Subsurface mesh does not contain material IDs");
    }
    if (subsurface_material_ids->size() % surface_material_ids->size() != 0)
    {
        OGS_FATAL("Could not determine the number of subsurface layers.");
    }
    int const number_of_layers =
        subsurface_material_ids->size() / surface_material_ids->size();

    ranges::copy(
        ranges::views::repeat_n(*surface_material_ids, number_of_layers) |
            ranges::views::join,
        subsurface_material_ids->begin());
}

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Creates a layered 3D OGS mesh from an existing 2D OGS mesh and a list "
        "of raster files representing subsurface layers. Supported raster "
        "formats are ArcGIS ascii rasters (*.asc), Surfer Grids (*.grd), or "
        "gridded XYZ rasters (*.xyz)."
        "Only input meshes consisting of line and triangle elements are "
        "currently supported as mapping of quads might result in invalid mesh "
        "elements.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2026, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::SwitchArg use_ascii_arg("", "ascii_output",
                                   "Write VTU output in ASCII format.");
    cmd.add(use_ascii_arg);

    double min_thickness(std::numeric_limits<double>::epsilon());
    TCLAP::ValueArg<double> min_thickness_arg(
        "t", "thickness",
        "The minimum thickness of a layer to be integrated at any given "
        "location, (min = 0)",
        false, min_thickness, "MIN_THICKNESS");
    cmd.add(min_thickness_arg);

    TCLAP::ValueArg<std::string> raster_path_arg(
        "r", "raster-list",
        "Input (.vtu). An ascii-file containing a list of input "
        "raster files, starting from"
        "top (DEM) to bottom",
        true, "", "INPUT_FILE_LIST");
    cmd.add(raster_path_arg);

    TCLAP::SwitchArg keep_materials_arg(
        "", "keep-surface-material-ids",
        "if the argument is present the materials defined in the surface mesh "
        "are used to set the material information for the subsurface cells",
        false);
    cmd.add(keep_materials_arg);

    TCLAP::ValueArg<std::string> mesh_out_arg(
        "o", "output-mesh-file",
        "Output (.vtu). The file name of the resulting 3D mesh", true, "",
        "OUTPUT_FILE");
    cmd.add(mesh_out_arg);

    TCLAP::ValueArg<std::string> mesh_arg(
        "i", "input-mesh-file",
        "Input (.vtu). The file name of the 2D input mesh", true, "",
        "INPUT_FILE");
    cmd.add(mesh_arg);
    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);

    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

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
    if (raster_paths.size() < 2)
    {
        ERR("At least two raster files needed to create 3D mesh.");
        return EXIT_FAILURE;
    }
    std::reverse(raster_paths.begin(), raster_paths.end());

    MeshToolsLib::MeshLayerMapper mapper;
    if (auto const rasters = FileIO::readRasters(raster_paths))
    {
        if (!mapper.createLayers(*sfc_mesh, *rasters, min_thickness))
        {
            return EXIT_FAILURE;
        }
    }
    else
    {
        ERR("Reading raster files.");
        return EXIT_FAILURE;
    }

    auto const result_mesh = mapper.getMesh("SubsurfaceMesh");
    if (result_mesh == nullptr)
    {
        ERR("Mapper returned empty result for 'SubsurfaceMesh'.");
        return EXIT_FAILURE;
    }

    if (keep_materials_arg.getValue())
    {
        assignSurfaceMaterialIDsToSubsurfaceLayers(sfc_mesh.get(),
                                                   result_mesh.get());
    }

    std::string output_name(mesh_out_arg.getValue());
    if (!BaseLib::hasFileExtension(".vtu", output_name))
    {
        output_name.append(".vtu");
    }

    INFO("Writing mesh '{:s}' ... ", output_name);
    auto const data_mode =
        use_ascii_arg.getValue() ? vtkXMLWriter::Ascii : vtkXMLWriter::Binary;

    MeshLib::IO::writeVtu(*result_mesh, output_name, data_mode);
    INFO("done.");

    return EXIT_SUCCESS;
}
