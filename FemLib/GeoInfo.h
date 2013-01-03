/**
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file GeoInfo.h
 *
 * Created on 2010-06-18 by Thomas Fischer
 */

#ifndef GEOINFO_H_
#define GEOINFO_H_

// GEO
#include "GeoObject.h"
#include "GeoType.h"

/**
 * \brief GeoInfo stores the type of the geometric entity and
 * the index within the vector the geometric entity is
 * managed. Possible geometric entities are documented in
 * GeoType.h.
 */
class GeoInfo
{
public:
	/**
	 * standard constructor. You need to set the attributes via
	 * setGeoType() and setGeoObj()!
	 */
	GeoInfo ();
	/**
	 * The constructor of a GeoInfo object initializes the
	 * attributes of the object.
	 * @param geo_type the type of the geometric entity.
	 * Possible geometric entities are documented in GeoType.h.
	 * @param geo_obj the pointer to an object of class GeoObject
	 */
	GeoInfo(GeoLib::GEOTYPE geo_type, const GeoLib::GeoObject* geo_obj = NULL);
	/**
	 * virtual destructor - destroys the object
	 */
	virtual ~GeoInfo();

	/**
	 * getter method for the geo type
	 * @sa enum GeoType
	 * @return the geo type
	 */
	GeoLib::GEOTYPE getGeomType () const;

	/**
	 * get the type as a string for log output
	 * @return
	 */
	std::string getGeomTypeAsString () const;

	/**
	 * getter for the pointer to the object
	 * @return
	 */
	const GeoLib::GeoObject* getGeoObj () const;

	/**
	 * setter for the geo type
	 * @sa enum GeoType
	 * @param geo_type type of the geometric entity
	 */
	void setGeoType (GeoLib::GEOTYPE geo_type);
	/**
	 * setter for the pointer to the GeoObject object
	 * @param geo_obj an instance of class GeoObject
	 */
	void setGeoObj (const GeoLib::GeoObject* geo_obj);

protected:
	/**
	 * type of the geometric entity. @sa enum GeoType
	 */
	GeoLib::GEOTYPE _geo_type;
	/**
	 * pointer to geometric object (GeoLib::Point, GeoLib::Polyline, GeoLib::Surface, ...)
	 */
	const GeoLib::GeoObject* _geo_obj;
};
#endif                                            /* GEOINFO_H_ */
