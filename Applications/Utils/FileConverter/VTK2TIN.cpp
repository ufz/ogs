// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <tclap/CmdLine.h>

#include <fstream>
#include <memory>
#include <string>

#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/TINInterface.h"
#include "GeoLib/Surface.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshToolsLib/convertMeshToGeo.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Converts VTK mesh into TIN file.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> mesh_in(
        "i", "mesh-input-file",
        "Input (.vtk). The name of the file containing the input mesh", true,
        "", "INPUT_FILE");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> mesh_out(
        "o", "TIN-output-file",
        "Output (.tin). The name of the file the TIN will be written to", true,
        "", "OUTPUT_FILE");
    cmd.add(mesh_out);
    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());
    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::VtuInterface::readVTUFile(mesh_in.getValue()));
    INFO("Mesh read: {:d} nodes, {:d} elements.", mesh->getNumberOfNodes(),
         mesh->getNumberOfElements());

    INFO("Converting the mesh to TIN");
    GeoLib::GEOObjects geo_objects;
    if (MeshToolsLib::convertMeshToGeo(*mesh, geo_objects))
    {
        INFO("Writing TIN into the file");
        GeoLib::IO::TINInterface::writeSurfaceAsTIN(
            *(*geo_objects.getSurfaceVec(mesh->getName()))[0],
            mesh_out.getValue());
    }

    return EXIT_SUCCESS;
}
