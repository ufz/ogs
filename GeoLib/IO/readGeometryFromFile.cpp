/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "readGeometryFromFile.h"

#include <vector>

#include "BaseLib/Error.h"
#include "BaseLib/FileTools.h"

#include "GeoLib/IO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "GeoLib/IO/Legacy/OGSIOVer4.h"

#include "GeoLib/GEOObjects.h"

namespace GeoLib
{
namespace IO
{
void
readGeometryFromFile(std::string const& fname, GeoLib::GEOObjects & geo_objs)
{
    if (BaseLib::getFileExtension(fname).compare("gml") == 0) {
        GeoLib::IO::BoostXmlGmlInterface xml(geo_objs);
        xml.readFile(fname);
    } else {
        std::vector<std::string> errors;
        std::string geo_name; // geo_name is output of the reading function
        GeoLib::IO::Legacy::readGLIFileV4(fname, geo_objs, geo_name, errors);
    }

    std::vector<std::string> geo_names;
    geo_objs.getGeometryNames(geo_names);
    if (geo_names.empty()) {
        OGS_FATAL("GEOObjects has no geometry name after reading the geometry file. "
                  "Something is wrong in the reading function.");
    }
}
}
}
