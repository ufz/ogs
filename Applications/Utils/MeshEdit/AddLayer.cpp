/*
 * \file
 * \brief Adds a layer to an existing mesh.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <memory>

#include "BaseLib/FileTools.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEditing/AddLayerToMesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Adds a layer to an existing mesh."
        "The documentation is available at "
        "https://docs.opengeosys.org/docs/tools/meshing/addlayer.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2022, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> mesh_arg(
        "i", "input-mesh-file", "the name of the file containing the mesh",
        true, "", "file name");
    cmd.add(mesh_arg);

    TCLAP::ValueArg<std::string> mesh_out_arg(
        "o", "output-mesh-file",
        "the name of the file the mesh should be written to (vtu format)", true,
        "", "file name");
    cmd.add(mesh_out_arg);

    TCLAP::ValueArg<double> layer_thickness_arg(
        "t", "layer-tickness", "the thickness of the new layer", false, 10,
        "floating point value");
    cmd.add(layer_thickness_arg);

    TCLAP::SwitchArg layer_position_arg(
        "", "add-layer-on-bottom",
        "Per default the layer is add on the top, if this argument is set the "
        "layer is add on the bottom.",
        true);
    cmd.add(layer_position_arg);

    TCLAP::SwitchArg copy_material_ids_arg(
        "", "copy-material-ids",
        "Copy the existing material distribution of the layer which is to be "
        "extended. If the switch isn't given a new material id will be "
        "created.",
        false);
    cmd.add(copy_material_ids_arg);

    cmd.parse(argc, argv);

    INFO("Reading mesh '{:s}' ... ", mesh_arg.getValue());
    auto subsfc_mesh = std::unique_ptr<MeshLib::Mesh>(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue()));
    if (!subsfc_mesh)
    {
        ERR("Error reading mesh '{:s}'.", mesh_arg.getValue());
        return EXIT_FAILURE;
    }
    INFO("done.");

    std::unique_ptr<MeshLib::Mesh> result(MeshLib::addLayerToMesh(
        *subsfc_mesh, layer_thickness_arg.getValue(), mesh_out_arg.getValue(),
        layer_position_arg.getValue(), copy_material_ids_arg.getValue()));
    if (!result)
    {
        ERR("Failure while adding layer.");
        return EXIT_FAILURE;
    }

    INFO("Writing mesh '{:s}' ... ", mesh_out_arg.getValue());
    MeshLib::IO::writeMeshToFile(*result, mesh_out_arg.getValue());
    INFO("done.");

    return EXIT_SUCCESS;
}
