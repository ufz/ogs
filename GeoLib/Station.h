/**
 * \file
 * \author Karsten Rink
 * \date   2010-07-01
 * \brief  Definition of the Station class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <string>

#include "Point.h"
#include "SensorData.h"

namespace GeoLib
{
/**
 * \brief A Station (observation site) is basically a Point with some additional
 * information.
 *
 * Additional information is largely optional (except for a name, but even this
 * may be empty). It may include a name, a stratigraphy (only for the derived
 * class StationBore), time series data from data loggers (as a
 * SensorData-object), etc.
 *
 * \sa StationBorehole, SensorData
 */
class Station : public Point
{
public:
    /**
     * \brief Constructor
     *
     * Constructor initialising a Station object
     * \param x The x-coordinate of the station.
     * \param y The y-coordinate of the station.
     * \param z The z-coordinate of the station.
     * \param name The name of the station.
     */
    explicit Station(double x = 0.0, double y = 0.0, double z = 0.0,
                     std::string name = "");

    explicit Station(Point const* coords, std::string name = "");

    /**
     * Constructor copies the source object
     * @param src the Station object that should be copied
     */
    Station(Station const& src);

    /// Returns the name of the station.
    std::string const& getName() const { return _name; }

    void setName(std::string const& name) { _name = name; }

    /// Creates a Station-object from information contained in a string
    /// (assuming the string has the right format)
    static Station* createStation(const std::string& line);

    /// Creates a new station object based on the given parameters.
    static Station* createStation(const std::string& name, double x, double y,
                                  double z);

    /// Returns the specific value for this station
    double getStationValue() const { return this->_station_value; }

    /// Allows to set a specific value for this station (e.g. for
    /// classification)
    void setStationValue(double station_value)
    {
        this->_station_value = station_value;
    }

    /// Allows to add sensor data from a CSV file to the observation site
    void addSensorDataFromCSV(const std::string& file_name)
    {
        _sensor_data.reset(new SensorData(file_name));
    }

    /// Returns all the sensor data for this observation site
    const SensorData* getSensorData() const { return _sensor_data.get(); }

private:
    std::string _name;
    double _station_value{0.0};
    std::unique_ptr<SensorData> _sensor_data{nullptr};
};

bool isStation(GeoLib::Point const* pnt);
}  // namespace GeoLib
