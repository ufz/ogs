/**
 * \file
 * \author Karsten Rink
 * \date   2010-01-25
 * \brief  Implementation of the SHPInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * @file SHPInterface.h
 * @date 25/01/2010
 * @author Karsten Rink
 */

#ifndef SHPINTERFACE_H
#define SHPINTERFACE_H

#include <string>

//ShapeLib includes
#include "shapefil.h"

// GeoLib
#include "GEOObjects.h"

/**
 * \brief Manages the import of ESRI shape files into GeoLib.
 */
class SHPInterface
{
public:
	/// Connection between ESRI type system for shape files and OGS GeoLib.
	enum OGSType
	{
		UNDEFINED   = 0,
		POINT       = 1,
		STATION     = 2,
		POLYLINE    = 3,
		POLYGON     = 4
	};

	/// Constructor
	SHPInterface(GeoLib::GEOObjects* geoObjects) : _geoObjects(geoObjects) {}

	/// Reads the header of the shape file.
	bool readSHPInfo(const std::string &filename, int &shapeType, int &numberOfEntities);

	/// Reads data from the shape file.
	void readSHPFile(const std::string &filename, OGSType choice, std::string listName);

private:
	/// Reads points into a vector of Point objects.
	void readPoints    (const SHPHandle &hSHP, int numberOfElements, std::string listName);

	/// Reads points into a vector of Point objects and marks them as Station.
	void readStations  (const SHPHandle &hSHP, int numberOfElements, std::string listName);

	/// Reads lines into a vector of Polyline objects.
	void readPolylines (const SHPHandle &hSHP, int numberOfElements, std::string listName);

	/// Reads lines into a vector of Polyline and Surface objects.
	void readPolygons  (const SHPHandle &hSHP, int numberOfElements, std::string listName);

	void adjustPolylines (std::vector<GeoLib::Polyline*>* lines,
	                      std::vector<std::size_t>  id_map);

	GeoLib::GEOObjects* _geoObjects;
};

#endif //SHPINTERFACE_H
