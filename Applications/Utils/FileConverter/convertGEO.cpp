/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <tclap/CmdLine.h>

#include <string>
#include <vector>

#include "Applications/FileIO/readGeometryFromFile.h"
#include "Applications/FileIO/writeGeometryToFile.h"
#include "BaseLib/MPI.h"
#include "GeoLib/GEOObjects.h"
#include "InfoLib/GitInfo.h"

int main(int argc, char* argv[])
{
    TCLAP::CmdLine cmd(
        "Converts OGS geometry file into another file format. "
        "Currently *.gml (OGS6 XML-based format) and *.gli (OGS5 format) "
        "formats are supported.\n\n"
        "OpenGeoSys-6 software, version " +
            GitInfoLib::GitInfo::ogs_version +
            ".\n"
            "Copyright (c) 2012-2025, OpenGeoSys Community "
            "(http://www.opengeosys.org)",
        ' ', GitInfoLib::GitInfo::ogs_version);
    TCLAP::ValueArg<std::string> argInputFileName(
        "i", "input-file", "the name of the geometry file to be converted",
        true, "", "file name");
    cmd.add(argInputFileName);
    TCLAP::ValueArg<std::string> argOutputFileName(
        "o", "output-file",
        "the name of the new geometry file whose file format is guessed from "
        "its file extension",
        true, "", "file name");
    cmd.add(argOutputFileName);

    TCLAP::ValueArg<std::string> gmsh_path_arg("g", "gmsh-path",
                                               "the path to the gmsh binary",
                                               false, "", "path as string");
    cmd.add(gmsh_path_arg);
    cmd.parse(argc, argv);

    BaseLib::MPI::Setup mpi_setup(argc, argv);

    GeoLib::GEOObjects geoObjects;
    FileIO::readGeometryFromFile(argInputFileName.getValue(), geoObjects,
                                 gmsh_path_arg.getValue());
    auto const geo_names = geoObjects.getGeometryNames();
    assert(geo_names.size() == 1);

    FileIO::writeGeometryToFile(geo_names[0], geoObjects,
                                argOutputFileName.getValue());

    return EXIT_SUCCESS;
}
