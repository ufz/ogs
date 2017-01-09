/*
 * \date 2015-04-14
 * \brief Adds a top layer to an existing mesh.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <memory>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "BaseLib/FileTools.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/AddLayerToMesh.h"

int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd(
        "Adds a top layer to an existing mesh"
        "The documentation is available at https://docs.opengeosys.org/docs/tools/meshing/addtoplayer ",
        ' ',
        "0.1");

    TCLAP::ValueArg<std::string> mesh_arg("i", "input-mesh-file",
        "the name of the file containing the mesh", true,
        "", "file name");
    cmd.add(mesh_arg);

    TCLAP::ValueArg<std::string> mesh_out_arg("o", "output-mesh-file",
        "the name of the file the mesh should be written to (vtu format)", true,
        "", "file name");
    cmd.add(mesh_out_arg);

    TCLAP::ValueArg<double> layer_thickness_arg("t", "layer-tickness",
        "the thickness of the new layer", false, 10, "floating point value");
    cmd.add(layer_thickness_arg);

    cmd.parse(argc, argv);

    INFO("Reading mesh \"%s\" ... ", mesh_arg.getValue().c_str());
    auto subsfc_mesh = std::unique_ptr<MeshLib::Mesh>(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));
    if (!subsfc_mesh) {
        ERR("Error reading mesh \"%s\".", mesh_arg.getValue().c_str());
        return EXIT_FAILURE;
    }
    INFO("done.");

    std::unique_ptr<MeshLib::Mesh> result(MeshLib::addTopLayerToMesh(
        *subsfc_mesh, layer_thickness_arg.getValue(), mesh_out_arg.getValue()));
    if (!result) {
        ERR("Failure while adding top layer.")
        return EXIT_FAILURE;
    }

    INFO("Writing mesh \"%s\" ... ", mesh_out_arg.getValue().c_str());
    MeshLib::IO::writeMeshToFile(*result, mesh_out_arg.getValue());
    INFO("done.");

    return EXIT_SUCCESS;
}
