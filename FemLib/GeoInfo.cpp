/*
 * GeoInfo.cpp
 *
 *  Created on: Jun 18, 2010
 *      Author: TF
 */

// STL
#include <limits>

// FEM
#include "GeoInfo.h"

GeoInfo::GeoInfo() :
	_geo_type(GEOLIB::INVALID), _geo_obj(NULL)
{}

GeoInfo::GeoInfo(GEOLIB::GEOTYPE geo_type, const GEOLIB::GeoObject* geo_obj) :
	_geo_type(geo_type), _geo_obj(geo_obj)
{}

GeoInfo::~GeoInfo()
{}

GEOLIB::GEOTYPE GeoInfo::getGeoType () const
{
	return _geo_type;
}

std::string GeoInfo::getGeoTypeAsString () const
{
	switch (_geo_type)
	{
	case GEOLIB::POINT:
		return "POINT";
	case GEOLIB::POLYLINE:
		return "POLYLINE";
	case GEOLIB::SURFACE:
		return "SURFACE";
	case GEOLIB::VOLUME:
		return "VOLUME";
	case GEOLIB::GEODOMAIN:
		return "DOMAIN";
	default:
		return "";
	}
}

void GeoInfo::setGeoType (GEOLIB::GEOTYPE geo_type)
{
	_geo_type = geo_type;
}

const GEOLIB::GeoObject* GeoInfo::getGeoObj () const
{
	return _geo_obj;
}

void GeoInfo::setGeoObj (const GEOLIB::GeoObject* geo_obj)
{
	_geo_obj = geo_obj;
}
