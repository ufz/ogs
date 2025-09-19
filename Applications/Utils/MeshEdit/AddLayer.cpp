/*
 * \file
 * \brief Adds a layer to an existing mesh.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include <memory>

#include "BaseLib/FileTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/MeshEditing/AddLayerToMesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Adds a layer to an existing mesh."
        "The documentation is available at "
        "https://www.opengeosys.org/docs/tools/meshing/addlayer."
        "\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(https://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> mesh_arg("i", "input-mesh-file",
                                          "Input (.vtu). The name of the file "
                                          "containing the mesh",
                                          true, "", "INPUT_FILE");
    cmd.add(mesh_arg);

    TCLAP::ValueArg<std::string> mesh_out_arg(
        "o", "output-mesh-file",
        "Output (.vtu). The name of the file the mesh should be written to",
        true, "", "OUTPUT_FILE");
    cmd.add(mesh_out_arg);

    TCLAP::ValueArg<double> layer_thickness_arg(
        "t", "layer-tickness",
        "the thickness of the new layer, "
        "(min = 0)",
        false, 10, "LAYER_THICKNESS");
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

    TCLAP::ValueArg<int> set_material_arg(
        "", "set-material-id", "the material id of the new layer, (min = 0)",
        false, 0, "MATERIAL_ID");
    cmd.add(set_material_arg);
    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    if (set_material_arg.isSet() && copy_material_ids_arg.isSet())
    {
        ERR("It is not possible to set both options '--copy-material-ids' and "
            "'--set-material-id'.");
        return EXIT_FAILURE;
    }

    INFO("Reading mesh '{:s}' ... ", mesh_arg.getValue());
    auto subsfc_mesh = std::unique_ptr<MeshLib::Mesh>(
        MeshLib::IO::readMeshFromFile(mesh_arg.getValue(), true));
    if (!subsfc_mesh)
    {
        ERR("Error reading mesh '{:s}'.", mesh_arg.getValue());
        return EXIT_FAILURE;
    }
    INFO("done.");

    std::optional<int> layer_id;
    if (set_material_arg.isSet())
    {
        layer_id = set_material_arg.getValue();
    }

    std::unique_ptr<MeshLib::Mesh> result(MeshToolsLib::addLayerToMesh(
        *subsfc_mesh, layer_thickness_arg.getValue(), mesh_out_arg.getValue(),
        layer_position_arg.getValue(), copy_material_ids_arg.getValue(),
        layer_id));
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
