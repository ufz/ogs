/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Applications/FileIO/GocadIO/GocadSGridReader.h"

#include <spdlog/spdlog.h>
#include <tclap/CmdLine.h>

#include <fstream>
#include <sstream>
#include <string>

#include "Applications/FileIO/GocadIO/GenerateFaceSetMeshes.h"
#include "BaseLib/FileTools.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshInformation.h"
#include "MeshLib/Node.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Reads a Gocad stratigraphic grid file (file ending sg) and writes the "
        "data in the vtk unstructured grid file format. The documentation is "
        "available at  "
        "https://www.opengeosys.org/docs/tools/meshing/gocadsgridreader/.\n\n "
        "OpenGeoSys-6 "
        "software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2021, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::SwitchArg face_set_arg("f", "generate-face-sets",
                                  "Generate face sets.");
    cmd.add(face_set_arg);

    TCLAP::ValueArg<std::string> mesh_output_arg(
        "o", "output-mesh", "vtk unstructured grid file name", true, "",
        "file_name.vtu");
    cmd.add(mesh_output_arg);

    TCLAP::ValueArg<std::string> sg_file_arg(
        "s", "sg", "Gocad stratigraphic grid file name", true, "",
        "file_name.sg");
    cmd.add(sg_file_arg);

    cmd.parse(argc, argv);

    // read the Gocad SGrid
    INFO("Start reading Gocad SGrid.");
    FileIO::Gocad::GocadSGridReader reader(sg_file_arg.getValue());
    INFO("End reading Gocad SGrid.");

    if (face_set_arg.getValue())
    {
        INFO("Generating a mesh for every face set.");
        FileIO::Gocad::generateFaceSets(
            reader, BaseLib::extractPath(sg_file_arg.getValue()));
    }

    std::unique_ptr<MeshLib::Mesh> mesh(reader.getMesh());
    MeshLib::MeshInformation::writePropertyVectorInformation(*mesh);

    INFO("Writing mesh to '{:s}'.", mesh_output_arg.getValue());
    MeshLib::IO::writeMeshToFile(*mesh, mesh_output_arg.getValue());

    return 0;
}
