/**
 * \file
 * \author Karsten Rink
 * \date   2010-07-01
 * \brief  Implementation of the Station class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>

// ThirdParty/logog
#include "logog/include/logog.hpp"

// BaseLib
#include "DateTools.h"
#include "StringTools.h"
// GeoLib
#include "Station.h"

namespace GeoLib
{
Station::Station(double x, double y, double z, std::string name) :
	Point (x,y,z), _name(name), _type(Station::StationType::STATION), _station_value(0.0), _sensor_data(NULL)
{}

Station::Station(Point* coords, std::string name) :
	Point (*coords), _name(name), _type(Station::StationType::STATION), _station_value(0.0), _sensor_data(NULL)
{}

Station::Station(Station const& src) :
	Point(src.getCoords()), _name(src._name), _type(src._type),
	_station_value(src._station_value), _sensor_data(NULL)
{}

Station::~Station()
{
	delete this->_sensor_data;
}

Station* Station::createStation(const std::string & line)
{
	std::list<std::string>::const_iterator it;
	Station* station = new Station();
	std::list<std::string> fields = BaseLib::splitString(line, '\t');

	if (fields.size() >= 3)
	{
		it = fields.begin();
		station->_name  = *it;
		(*station)[0]     = strtod((BaseLib::replaceString(",", ".", *(++it))).c_str(), NULL);
		(*station)[1]     = strtod((BaseLib::replaceString(",", ".", *(++it))).c_str(), NULL);
		if (++it != fields.end())
			(*station)[2] = strtod((BaseLib::replaceString(",", ".", *it)).c_str(), NULL);
	}
	else
	{
		INFO("Station::createStation() - Unexpected file format.");
		delete station;
		return NULL;
	}
	return station;
}

Station* Station::createStation(const std::string &name, double x, double y, double z)
{
	Station* station = new Station();
	station->_name = name;
	(*station)[0] = x;
	(*station)[1] = y;
	(*station)[2] = z;
	return station;
}

} // namespace

