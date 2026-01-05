// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <tclap/CmdLine.h>

#include "Applications/FileIO/SHPInterface.h"
#include "BaseLib/Logging.h"
#include "BaseLib/MPI.h"
#include "BaseLib/TCLAPArguments.h"
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
            "Copyright (c) 2012-2026, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);

    TCLAP::ValueArg<std::string> output_arg("o", "output-file",
                                            "Output (.shp). Esri Shapefile",
                                            true, "", "OUTPUT_FILE");
    cmd.add(output_arg);

    TCLAP::ValueArg<std::string> input_arg("i", "input-file",
                                           "Input (.vtu | .msh). OGS mesh file",
                                           true, "", "INPUT_FILE");
    cmd.add(input_arg);

    auto log_level_arg = BaseLib::makeLogLevelArg();
    cmd.add(log_level_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);
    BaseLib::initOGSLogger(log_level_arg.getValue());

    std::string const file_name(input_arg.getValue());
    std::unique_ptr<MeshLib::Mesh> const mesh(
        MeshLib::IO::readMeshFromFile(file_name));
    if (FileIO::SHPInterface::write2dMeshToSHP(output_arg.getValue(), *mesh))
    {
        return EXIT_SUCCESS;
    }
    return EXIT_FAILURE;
}
