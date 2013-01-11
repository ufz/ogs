/**
 * \file
 * \author Thomas Fischer
 * \date   2010-06-18
 * \brief  Implementation of the GeoInfo class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// STL
#include <limits>

// FEM
#include "GeoInfo.h"

GeoInfo::GeoInfo() :
	_geo_type(GeoLib::INVALID), _geo_obj(NULL)
{}

GeoInfo::GeoInfo(GeoLib::GEOTYPE geo_type, const GeoLib::GeoObject* geo_obj) :
	_geo_type(geo_type), _geo_obj(geo_obj)
{}

GeoInfo::~GeoInfo()
{}

GeoLib::GEOTYPE GeoInfo::getGeomType () const
{
	return _geo_type;
}

std::string GeoInfo::getGeomTypeAsString () const
{
	switch (_geo_type)
	{
	case GeoLib::POINT:
		return "POINT";
	case GeoLib::POLYLINE:
		return "POLYLINE";
	case GeoLib::SURFACE:
		return "SURFACE";
	case GeoLib::VOLUME:
		return "VOLUME";
	case GeoLib::GEODOMAIN:
		return "DOMAIN";
	default:
		return "";
	}
}

void GeoInfo::setGeoType (GeoLib::GEOTYPE geo_type)
{
	_geo_type = geo_type;
}

const GeoLib::GeoObject* GeoInfo::getGeoObj () const
{
	return _geo_obj;
}

void GeoInfo::setGeoObj (const GeoLib::GeoObject* geo_obj)
{
	_geo_obj = geo_obj;
}
