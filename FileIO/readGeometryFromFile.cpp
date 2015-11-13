/**
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "readGeometryFromFile.h"

#include <vector>

#include "BaseLib/FileTools.h"

#include "FileIO/XmlIO/Boost/BoostXmlGmlInterface.h"
#include "FileIO/Legacy/OGSIOVer4.h"

#include "GeoLib/GEOObjects.h"

namespace FileIO
{
void
readGeometryFromFile(std::string const& fname, GeoLib::GEOObjects & geo_objs)
{
	if (BaseLib::getFileExtension(fname).compare("gml") == 0) {
		FileIO::BoostXmlGmlInterface xml(geo_objs);
		xml.readFile(fname);
	} else {
		std::vector<std::string> errors;
		std::string geo_name(BaseLib::extractBaseNameWithoutExtension(fname));
		FileIO::Legacy::readGLIFileV4(fname, &geo_objs, geo_name, errors);
	}
}
}
