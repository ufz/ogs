/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <string>
#include <vector>

#include <tclap/CmdLine.h>

#include "Applications/ApplicationsLib/LogogSetup.h"

#include "GeoLib/IO/readGeometryFromFile.h"
#include "GeoLib/IO/writeGeometryToFile.h"
#include "GeoLib/GEOObjects.h"


int main (int argc, char* argv[])
{
    ApplicationsLib::LogogSetup logog_setup;

    TCLAP::CmdLine cmd("Converts OGS geometric file into another file format.", ' ', "0.1");
    TCLAP::ValueArg<std::string> argInputFileName("i", "input-file",
                                         "the name of the gli file to be converted", true,
                                         "", "file name");
    cmd.add(argInputFileName);
    TCLAP::ValueArg<std::string> argOutputFileName("o", "output-file",
                                          "the name of the new gml file", true,
                                          "", "file name");
    cmd.add(argOutputFileName);
    cmd.parse(argc, argv);

    GeoLib::GEOObjects geoObjects;
    GeoLib::IO::readGeometryFromFile(argInputFileName.getValue(), geoObjects);
    std::vector<std::string> geo_names;
    geoObjects.getGeometryNames(geo_names);
    assert(geo_names.size() == 1);

    GeoLib::IO::writeGeometryToFile(geo_names[0], geoObjects, argOutputFileName.getValue());

    return EXIT_SUCCESS;
}
