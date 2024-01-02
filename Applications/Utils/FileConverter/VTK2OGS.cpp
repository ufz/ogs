/**
 * \file VTK2OGS.cpp
 * \author Norihiro Watanabe
 * \date Aug 07, 2013
 * \brief Converts VTK mesh into OGS mesh.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

// STL
#include <tclap/CmdLine.h>

#ifdef USE_PETSC
#include <mpi.h>
#endif

#include <string>

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
            "Copyright (c) 2012-2024, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> mesh_in(
        "i", "mesh-input-file",
        "the name of the file containing the input mesh", true, "",
        "file name of input mesh");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> mesh_out(
        "o", "mesh-output-file",
        "the name of the file the mesh will be written to", true, "",
        "file name of output mesh");
    cmd.add(mesh_out);
    cmd.parse(argc, argv);

#ifdef USE_PETSC
    MPI_Init(&argc, &argv);
#endif

    MeshLib::Mesh* mesh(
        MeshLib::IO::VtuInterface::readVTUFile(mesh_in.getValue()));
    INFO("Mesh read: {:d} nodes, {:d} elements.", mesh->getNumberOfNodes(),
         mesh->getNumberOfElements());

    MeshLib::IO::Legacy::MeshIO meshIO;
    meshIO.setMesh(mesh);
    BaseLib::IO::writeStringToFile(meshIO.writeToString(), mesh_out.getValue());

#ifdef USE_PETSC
    MPI_Finalize();
#endif
    return EXIT_SUCCESS;
}
