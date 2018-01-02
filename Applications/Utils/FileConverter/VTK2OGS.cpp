/**
 * @file VTK2OGS.cpp
 * @author Norihiro Watanabe
 * @date Aug 07, 2013
 * @brief Converts VTK mesh into OGS mesh.
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// STL
#include <string>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "MeshLib/IO/VtkIO/VtuInterface.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/IO/Legacy/MeshIO.h"

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Converts VTK mesh into OGS mesh.", ' ', "0.1");
    TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
                                         "the name of the file containing the input mesh", true,
                                         "", "file name of input mesh");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
                                          "the name of the file the mesh will be written to", true,
                                          "", "file name of output mesh");
    cmd.add(mesh_out);
    cmd.parse(argc, argv);

    MeshLib::Mesh* mesh (MeshLib::IO::VtuInterface::readVTUFile(mesh_in.getValue()));
    INFO("Mesh read: %d nodes, %d elements.", mesh->getNumberOfNodes(), mesh->getNumberOfElements());

    MeshLib::IO::Legacy::MeshIO meshIO;
    meshIO.setMesh(mesh);
    meshIO.writeToFile(mesh_out.getValue());

    return EXIT_SUCCESS;
}
