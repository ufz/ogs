/**
 * \file
 * \author Karsten Rink
 * \date   2010-07-01
 * \brief  Definition of the Station class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef GEO_STATION_H
#define GEO_STATION_H

#include <list>
#include <map>
#include <string>
#include <vector>

#include "Point.h"
#include "Polyline.h"
#include "PropertyBounds.h"
#include "SensorData.h"

namespace GeoLib
{
/**
 * \ingroup GeoLib
 *
 * \brief A Station (observation site) is basically a Point with some additional information.
 *
 * Additional information is largely optional (except for a name, but even this may be empty).
 * It may include a name, a stratigraphy (only for the derived class StationBore),
 * time series data from data loggers (as a SensorData-object), etc.
 *
 * Notes concerning the property-system used in this class:
 * Variables of Station and derived classes can be defined to be "properties" of this class (this is entirely optional!).
 * Certain functions in the GUI allow you to modify aspects of the visualisation based on these
 * properties (e.g. filtering operations such as "display only boreholes drilled after 1990 with a
 * depth between 400-800m").
 * To make use of this functionality you need to define properties using the "Station::addProperty()"-method.
 * Parameters for this function include the name of the property as well as a method to read and write the
 * value of the associated variable (see documentation for "addProperty" for details). Furthermore, these read
 * and write-functions need to be actually implemented as static functions to avoid casting problems with the
 * function pointers used to dynamically connect the GUI functionality to the variables defined within the
 * station-classes. Please refer to the documentation of the properties defined below for details.
 *
 * \sa StationBorehole, SensorData
 */
class Station : public Point
{
protected:

	//typedef double (Station::*getFct)();
	//typedef void (Station::*setFct)(double);

	/**
	 * \brief Container for station-properties.
	 * Each property consists of a name, a get- and a set-function.
	 * Please refer to Station::addProperty for details.
	 */
	struct STNProperty
	{
		std::string name;
		double (* get)(void*);
		void (* set)(void*, double);
	};

public:
	/// Signals if the object is a "simple" Station or a Borehole (i.e. containing borehole-specific information).
	enum class StationType
	{
		INVALID = 0,
		STATION,
		BOREHOLE
	};

	/**
	 * \brief Constructor
	 *
	 * Constructor initialising a Station object
	 * \param x The x-coordinate of the station.
	 * \param y The y-coordinate of the station.
	 * \param z The z-coordinate of the station.
	 * \param name The name of the station.
	 */
	Station(double x = 0.0,
	        double y = 0.0,
	        double z = 0.0,
	        std::string name = "");

	Station(Point* coords, std::string name = "");

	/**
	 * Constructor copies the source object
	 * @param src the Station object that should be copied
	 */
	Station(Station const& src);

	virtual ~Station();

	/**
	 * \brief Defines a property for this class.
	 *
	 * Variables in Station and its derived classes can be defined to be properties of this station.
	 * This definition consists of a name for the property as well as a function pointer to a method to
	 * read and to write this variable.
	 * Due to inheritance these function pointers only work correctly if the read and write functions are
	 * implemented as static functions, i.e. both the read and the write functin get a void-pointer to the
	 * actual station-object. This pointer is then casted to the correct class and the respective value can
	 * then be read or written dynamically from that object. It is highly recommended to define both the
	 * read and write function as protected because it does not actually make sense for these functions to be
	 * static except in the context of function pointers. Please refer to the examples below, i.e. the getX
	 * and setX methods.
	 * \param pname The name of the property.
	 * \param get A function pointer to a static read function for the variable referred to by pname
	 * \param set A function pointer to a static write function for the variable referred to by pname
	 */
	void addProperty(const std::string &pname, double (* get)(void*), void (* set)(void*, double));

	/// Returns a map containing all the properties of that station type.
	const std::map<std::string, double> getProperties();

	/// Determines if the station's parameters are within the the bounds of the current selection (see property system for details)
	bool inSelection(const std::vector<PropertyBounds> &bounds);

	/// Returns true if all properties of this stations are within the boundaries given by \param bounds and false otherwise
	bool inSelection(std::map<std::string, double> properties) const;

	/// Returns the name of the station.
	std::string const& getName() const { return _name; }

	/// Returns the GeoSys-station-type for the station.
	StationType type() const { return _type; }

	/// Creates a Station-object from information contained in a string (assuming the string has the right format)
	static Station* createStation(const std::string &line);

	/// Creates a new station object based on the given parameters.
	static Station* createStation(const std::string &name, double x, double y, double z);

	/// Returns the specific value for this station
	double getStationValue() { return this->_station_value; }

	/// Allows to set a specific value for this station (e.g. for classification)
	void setStationValue(double station_value) { this->_station_value = station_value; }

	/// Allows to add sensor data from a CSV file to the observation site
	void addSensorDataFromCSV(const std::string &file_name) { this->_sensor_data = new SensorData(file_name); }

	/// Returns all the sensor data for this observation site
	const SensorData* getSensorData() { return this->_sensor_data; }

protected:
	/**
	 * \brief Returns the x-coordinate of this station. See the detailed documentation for getX() concerning the syntax.
	 *
	 * Returns the x-coordinate of this station.
	 * This function belongs to the property system of Station and return the value for property "x"
	 * (i.e. the x-coordinate of the station). It is implemented as a static method to avoid casting issues
	 * related to the function pointer associated with this function. Therefore, this function needs to be
	 * called "getX((void*)this);". It is highly recommended to define this function as protected because it
	 * does not actually make sense for these functions to be static except in the context of function pointers.
	 * \param stnObject A pointer to the station object for which the x-coordinate should be returned, usually (void*)this will work fine.
	 * \return The x-coordinate for this station.
	 */
	static double getX(void* stnObject) { Station* stn = (Station*)stnObject;
		                              return (*stn)[0]; }
	/// Returns the y-coordinate of this station. See the detailed documentation for getX concerning the syntax.
	static double getY(void* stnObject) { Station* stn = (Station*)stnObject;
		                              return (*stn)[1]; }
	/// Returns the z-coordinate of this station. See the detailed documentation for getX concerning the syntax.
	static double getZ(void* stnObject) { Station* stn = (Station*)stnObject;
		                              return (*stn)[2]; }
	/// Sets the x-coordinate for this station. See the detailed documentation for getX concerning the syntax.
	static void setX(void* stnObject, double val) { Station* stn = (Station*)stnObject;
		                                        (*stn)[0] = val; }
	/// Sets the y-coordinate for this station. See the detailed documentation for getX concerning the syntax.
	static void setY(void* stnObject, double val) { Station* stn = (Station*)stnObject;
		                                        (*stn)[1] = val; }
	/// Sets the z-coordinate for this station. See the detailed documentation for getX concerning the syntax.
	static void setZ(void* stnObject, double val) { Station* stn = (Station*)stnObject;
		                                        (*stn)[2] = val; }

	std::string _name;
	StationType _type; // GeoSys Station Type
	std::vector<STNProperty> _properties;

private:
	double _station_value;
	SensorData* _sensor_data;

};

} // namespace

#endif // GEO_STATION_H
