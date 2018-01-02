/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <string>
#include <vector>

#include <tclap/CmdLine.h>

#include "BaseLib/BuildInfo.h"

#include "Applications/ApplicationsLib/LogogSetup.h"
#include "Applications/FileIO/readGeometryFromFile.h"
#include "Applications/FileIO/writeGeometryToFile.h"
#include "GeoLib/GEOObjects.h"


int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Converts OGS geometry file into another file format. "
                       "Currently *.gml (OGS6 XML-based format) and *.gli (OGS5 format) formats are supported.",
                       ' ', BaseLib::BuildInfo::git_describe);
    TCLAP::ValueArg<std::string> argInputFileName("i", "input-file",
                                         "the name of the geometry file to be converted", true,
                                         "", "file name");
    cmd.add(argInputFileName);
    TCLAP::ValueArg<std::string> argOutputFileName("o", "output-file",
                                          "the name of the new geometry file whose file format is guessed from its file extension", true,
                                          "", "file name");
    cmd.add(argOutputFileName);
    cmd.parse(argc, argv);

    GeoLib::GEOObjects geoObjects;
    FileIO::readGeometryFromFile(argInputFileName.getValue(), geoObjects);
    std::vector<std::string> geo_names;
    geoObjects.getGeometryNames(geo_names);
    assert(geo_names.size() == 1);

    FileIO::writeGeometryToFile(geo_names[0], geoObjects, argOutputFileName.getValue());

    return EXIT_SUCCESS;
}
