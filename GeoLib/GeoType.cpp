/**
 * \file
 * \author Thomas Fischer
 * \date   2010-12-01
 * \brief  Implementation of GEOTYPE enumeration helper functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GeoType.h"

namespace GeoLib {

GEOTYPE convertGeoType (const std::string& geo_type_str)
{
	if (geo_type_str.compare ("POINT") == 0) return GEOTYPE::POINT;
	if (geo_type_str.compare ("POLYLINE") == 0) return GEOTYPE::POLYLINE;
	if (geo_type_str.compare ("SURFACE") == 0) return GEOTYPE::SURFACE;
	if (geo_type_str.compare ("VOLUME") == 0) return GEOTYPE::VOLUME;
	if (geo_type_str.compare ("GEODOMAIN") == 0) return GEOTYPE::GEODOMAIN;
	return GEOTYPE::INVALID;
}

std::string convertGeoTypeToString (GEOTYPE geo_type)
{
	if (geo_type == GEOTYPE::POINT) return "POINT";
	if (geo_type == GEOTYPE::POLYLINE) return "POLYLINE";
	if (geo_type == GEOTYPE::SURFACE) return "SURFACE";
	if (geo_type == GEOTYPE::VOLUME) return "VOLUME";
	if (geo_type == GEOTYPE::GEODOMAIN) return "GEODOMAIN";
	return "INVALID";
}

} // end namespace GeoLib
