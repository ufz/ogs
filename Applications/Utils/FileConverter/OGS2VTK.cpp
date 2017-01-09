/**
 * @file OGS2VTK.cpp
 * @author Thomas Fischer
 * @date Jan 24, 2013
 * @brief Converts OGS mesh into VTK mesh.
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include <string>
#include <memory>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Converts OGS mesh into VTK mesh.", ' ', "0.1");
    TCLAP::ValueArg<std::string> mesh_in("i", "mesh-input-file",
                                         "the name of the file containing the input mesh", true,
                                         "", "file name of input mesh");
    cmd.add(mesh_in);
    TCLAP::ValueArg<std::string> mesh_out("o", "mesh-output-file",
                                          "the name of the file the mesh will be written to", true,
                                          "", "file name of output mesh");
    cmd.add(mesh_out);
    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh const> mesh(
        MeshLib::IO::readMeshFromFile(mesh_in.getValue()));
    INFO("Mesh read: %d nodes, %d elements.", mesh->getNumberOfNodes(), mesh->getNumberOfElements());

    MeshLib::IO::VtuInterface vtu(mesh.get());
    vtu.writeToFile(mesh_out.getValue());

    return EXIT_SUCCESS;
}
