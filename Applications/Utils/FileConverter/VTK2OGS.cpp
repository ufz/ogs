/**
 * \file VTK2OGS.cpp
 * \author Norihiro Watanabe
 * \date Aug 07, 2013
 * \brief Converts VTK mesh into OGS mesh.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// STL
#include <tclap/CmdLine.h>

#include <string>

#include "BaseLib/MPI.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/Legacy/MeshIO.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/Mesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Converts VTK mesh into OGS mesh.\n\n"
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
        "o", "mesh-output-file",
        "Output (.msh). The name of the file the mesh will be written to", true,
        "", "OUTPUT_FILE");
    cmd.add(mesh_out);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);

    MeshLib::Mesh* mesh(
        MeshLib::IO::VtuInterface::readVTUFile(mesh_in.getValue()));
    INFO("Mesh read: {:d} nodes, {:d} elements.", mesh->getNumberOfNodes(),
         mesh->getNumberOfElements());

    MeshLib::IO::Legacy::MeshIO meshIO;
    meshIO.setMesh(mesh);
    BaseLib::IO::writeStringToFile(meshIO.writeToString(), mesh_out.getValue());

    return EXIT_SUCCESS;
}
