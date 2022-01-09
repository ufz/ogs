/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <tclap/CmdLine.h>

#include "Applications/FileIO/SHPInterface.h"
#include "InfoLib/GitInfo.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Converts 2D mesh file into shapfile such that each element is "
        "represented by a polygon. Cell attributes are transferred onto shape "
        "polygons while point attributes are ignored.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2022, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> output_arg("o", "output-file",
                                            "Esri Shapefile (*.shp)", true, "",
                                            "output_file.shp");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> input_arg("i", "input-file",
                                           "OGS mesh file (*.vtu, *.msh)", true,
                                           "", "input_file.vtu");
    cmd.add(input_arg);

    cmd.parse(argc, argv);

    std::string const file_name(input_arg.getValue());
    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshLib::IO::readMeshFromFile(file_name));
    if (FileIO::SHPInterface::write2dMeshToSHP(output_arg.getValue(), *mesh))
    {
        return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
}
