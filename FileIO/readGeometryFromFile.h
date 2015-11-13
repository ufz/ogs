/**
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef READGEOMETRYFROMFILE_H
#define READGEOMETRYFROMFILE_H

#include <string>

namespace GeoLib
{
	class GEOObjects;
}

namespace FileIO
{
	void
	readGeometryFromFile(std::string const& fname, GeoLib::GEOObjects & geo_objs);
}

#endif // READGEOMETRYFROMFILE_H
