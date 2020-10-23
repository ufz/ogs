/**
 * \file
 * \author Thomas Fischer
 * \date Jan 24, 2013
 * \brief Converts OGS mesh into VTK mesh.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <string>
#include <memory>

#include <tclap/CmdLine.h>

#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"

int main (int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Converts OGS mesh into VTK mesh.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2020, OpenGeoSys Community "
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
    TCLAP::SwitchArg use_ascii_arg(
        "", "ascii_output",
        "Write VTU output in ASCII format. Due to possible rounding the ascii "
        "output could result in lower accuracy.");
    cmd.add(use_ascii_arg);
    cmd.parse(argc, argv);

    std::unique_ptr<MeshLib::Mesh const> mesh(
        MeshLib::IO::readMeshFromFile(mesh_in.getValue()));
    if (!mesh)
    {
        return EXIT_FAILURE;
    }
    INFO("Mesh read: {:d} nodes, {:d} elements.", mesh->getNumberOfNodes(),
         mesh->getNumberOfElements());

    auto const data_mode =
        use_ascii_arg.getValue() ? vtkXMLWriter::Ascii : vtkXMLWriter::Binary;

    MeshLib::IO::writeVtu(*mesh, mesh_out.getValue(), data_mode);

    return EXIT_SUCCESS;
}
