/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 *
 * \file GeoType.cpp
 *
 * Created on 2010-12-01 by Thomas Fischer
 */

#include "GeoType.h"

namespace GeoLib {

GEOTYPE convertGeoType (const std::string& geo_type_str)
{
	if (geo_type_str.compare ("POINT") == 0) return POINT;
	if (geo_type_str.compare ("POLYLINE") == 0) return POLYLINE;
	if (geo_type_str.compare ("SURFACE") == 0) return SURFACE;
	if (geo_type_str.compare ("VOLUME") == 0) return VOLUME;
	if (geo_type_str.compare ("GEODOMAIN") == 0) return GEODOMAIN;
	return INVALID;
}

std::string convertGeoTypeToString (GEOTYPE geo_type)
{
	if (geo_type == POINT) return "POINT";
	if (geo_type == POLYLINE) return "POLYLINE";
	if (geo_type == SURFACE) return "SURFACE";
	if (geo_type == VOLUME) return "VOLUME";
	if (geo_type == GEODOMAIN) return "GEODOMAIN";
	return "INVALID";
}

} // end namespace GeoLib
