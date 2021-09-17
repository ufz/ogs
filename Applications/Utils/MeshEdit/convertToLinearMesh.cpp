/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshLib/MeshEditing/ConvertToLinearMesh.h"

#include <tclap/CmdLine.h>

#include <memory>
#include <string>

#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Convert a non-linear mesh to a linear mesh.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> input_arg(
        "i", "input-mesh-file", "input mesh file", true, "", "string");
    cmd.add(input_arg);
    TCLAP::ValueArg<std::string> output_arg(
        "o", "output-mesh-file", "output mesh file", true, "", "string");
    cmd.add(output_arg);

    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh> mesh(
        MeshLib::IO::readMeshFromFile(input_arg.getValue()));
    if (!mesh)
    {
        return EXIT_FAILURE;
    }
    if (!mesh->hasNonlinearElement())
    {
        ERR("The input mesh is linear. Exit.");
        return EXIT_FAILURE;
    }

    INFO("Converting to a linear order mesh");
    std::unique_ptr<MeshLib::Mesh> new_mesh(
        MeshLib::convertToLinearMesh(*mesh, mesh->getName() + "_linear"));

    INFO("Save the new mesh into a file");
    MeshLib::IO::writeMeshToFile(*new_mesh, output_arg.getValue());

    return EXIT_SUCCESS;
}
