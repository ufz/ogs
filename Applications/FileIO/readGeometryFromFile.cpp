/**
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "Legacy/OGSIOVer4.h"

#include "GeoLib/GEOObjects.h"

namespace FileIO
{
void readGeometryFromFile(std::string const& fname,
                          GeoLib::GEOObjects& geo_objs,
                          std::string const& gmsh_path)
{
    if (BaseLib::getFileExtension(fname) == ".gml")
    {
        GeoLib::IO::BoostXmlGmlInterface xml(geo_objs);
        xml.readFile(fname);
    } else {
        std::vector<std::string> errors;
        std::string geo_name; // geo_name is output of the reading function
        FileIO::Legacy::readGLIFileV4(
            fname, geo_objs, geo_name, errors, gmsh_path);
    }

    if (geo_objs.getGeometryNames().empty())
    {
        OGS_FATAL(
            "GEOObjects has no geometry name after reading the geometry file. "
            "Something is wrong in the reading function.");
    }
}
}  // namespace FileIO
