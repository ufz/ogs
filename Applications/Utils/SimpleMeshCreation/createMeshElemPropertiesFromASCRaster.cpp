/**
 * \file
 * \brief  Implementation of the createMeshElemPropertiesFromASCRaster tool.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <tclap/CmdLine.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <numeric>

#include "BaseLib/FileTools.h"
#include "BaseLib/MPI.h"
#include "BaseLib/quicksort.h"
#include "GeoLib/IO/AsciiRasterInterface.h"
#include "GeoLib/Raster.h"
#include "InfoLib/GitInfo.h"
#include "MathLib/MathTools.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/VtkIO/VtkMeshConverter.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshToolsLib/MeshEditing/Mesh2MeshPropertyInterpolation.h"
#include "MeshToolsLib/MeshGenerators/RasterToMesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Generates properties for mesh elements of an input mesh deploying a "
        "ASC raster file.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> out_mesh_arg(
        "o",
        "out-mesh",
        "the mesh is stored to a file of this name",
        false,
        "",
        "filename for mesh output");
    cmd.add(out_mesh_arg);

    TCLAP::SwitchArg refinement_raster_output_arg(
        "", "output-refined-raster",
        "Write refined raster to an additional, new ASC file.");
    cmd.add(refinement_raster_output_arg);

    TCLAP::ValueArg<unsigned> refinement_arg(
        "r",
        "refine",
        "refinement factor that raises the resolution of the raster data",
        false,
        1,
        "factor (default = 1)");
    cmd.add(refinement_arg);

    TCLAP::ValueArg<std::string> raster_arg("",
                                            "raster-file",
                                            "the name of the ASC raster file",
                                            true,
                                            "",
                                            "file name");
    cmd.add(raster_arg);

    TCLAP::ValueArg<std::string> property_arg(
        "p",
        "property-name",
        "the name of the property the values are stored for",
        true,
        "",
        "property name as string");
    cmd.add(property_arg);

    TCLAP::ValueArg<std::string> mesh_arg(
        "m", "mesh", "the mesh is read from this file", true, "", "file name");
    cmd.add(mesh_arg);

    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);

    // read mesh
    std::unique_ptr<MeshLib::Mesh> dest_mesh(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));

    // read raster and if required manipulate it
    auto raster = std::unique_ptr<GeoLib::Raster>(
        FileIO::AsciiRasterInterface::getRasterFromASCFile(
            raster_arg.getValue()));
    GeoLib::RasterHeader header(raster->getHeader());
    if (refinement_arg.getValue() > 1)
    {
        raster->refineRaster(refinement_arg.getValue());
        if (refinement_raster_output_arg.getValue())
        {
            // write new asc file
            std::string new_raster_fname(
                BaseLib::dropFileExtension(raster_arg.getValue()));
            new_raster_fname += "-" + std::to_string(header.n_rows) + "x" +
                                std::to_string(header.n_cols) + ".asc";
            FileIO::AsciiRasterInterface::writeRasterAsASC(*raster,
                                                           new_raster_fname);
        }
    }

    std::unique_ptr<MeshLib::Mesh> src_mesh(
        MeshToolsLib::RasterToMesh::convert(*raster,
                                            MeshLib::MeshElemType::QUAD,
                                            MeshLib::UseIntensityAs::DATAVECTOR,
                                            property_arg.getValue()));

    // do the interpolation
    MeshToolsLib::Mesh2MeshPropertyInterpolation mesh_interpolation(
        *src_mesh, property_arg.getValue());
    mesh_interpolation.setPropertiesForMesh(*dest_mesh);

    if (!out_mesh_arg.getValue().empty())
    {
        MeshLib::IO::writeMeshToFile(*dest_mesh, out_mesh_arg.getValue());
    }

    return EXIT_SUCCESS;
}
