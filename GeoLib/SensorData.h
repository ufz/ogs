/**
 * \file
 * \author Karsten Rink
 * \date   2012-08-01
 * \brief  Definition of the SensorData class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <cstddef>
#include <string>
#include <vector>

/**
 * Possible types of input data for time series sensor data.
 * Implementation as Enum for specific implementations later on.
 *
 * \sa SensorData
 */
enum class SensorDataType
{
    OTHER = 0,
    PRECIPITATION,
    EVAPORATION,
    TEMPERATURE
    // please expand if necessary
};

/**
 * Possible types of time specification.
 * In addition to the usual units we added 'DATE' for specification of dates
 * in the format 'dd.mm.yyyy' as well as 'DATETIME' in the format
 * 'dd.mm.yyyy.hh.mm.ss'.
 */
enum class TimeStepType
{
    NONE = 0,
    SECONDS,
    MINUTES,
    DAYS,
    WEEKS,
    MONTHS,
    YEARS,
    DATE,    // time series is given as a vector of dates
    DATETIME // time series is given as a vector of date + time
};

/**
 * \brief A container for sensor data at an observation site.
 * The class stores a number of time series and has been created for use in Station-objects.
 *
 * \sa Station
 */
class SensorData final
{
public:
    /// Constructor using file name (automatically reads the file and fills all data structures)
    explicit SensorData(const std::string& file_name);

    /// Constructor using a time step vector valid for all time series that will be added later
    explicit SensorData(std::vector<std::size_t> time_steps);

    /// Constructor using time step bounds for all time series that will be added later
    SensorData(std::size_t first_timestep, std::size_t last_timestep, std::size_t step_size);

    ~SensorData();

    /// Adds a time series that needs to conform to the time step vector specified in the constructor.
    /// Optionally a unit for the time series can be given.
    /// The name is converted to SensorDataType enum.
    void addTimeSeries( const std::string &data_name, std::vector<float> *data, const std::string &data_unit_string = "" );

    /// Adds a time series that needs to conform to the time step vector specified in the constructor.
    /// Optionally a unit for the time series can be given.
    void addTimeSeries( SensorDataType data_name, std::vector<float> *data, const std::string &data_unit_string = "" );

    /// Returns the time series with the given name
    const std::vector<float>* getTimeSeries(SensorDataType time_series_name) const;

    /// Returns all time series names contained in this container
    const std::vector<SensorDataType>& getTimeSeriesNames() const { return vec_names_; }

    /// Returns the time step vector (if it exists)
    const std::vector<std::size_t>& getTimeSteps() const { return time_steps_; }

    /// Returns the first time step
    std::size_t getStartTime() const { return start_; }

    /// Returns the last time step
    std::size_t getEndTime() const { return end_; }

    /// Returns the interval between time steps (Returns "0" if a vector is given!)
    std::size_t getStepSize() const { return step_size_; }

    /// Allows to set a unit for the time steps
    void setTimeUnit(TimeStepType t) { time_unit_ = t; }

    /// Returns the unit the time steps
    TimeStepType getTimeUnit() const { return time_unit_; }

    /// Converts Sensor Data Types to Strings
    static std::string convertSensorDataType2String(SensorDataType t);

    /// Converts Strings to Sensor Data Types
    static SensorDataType convertString2SensorDataType(const std::string &s);

private:
    /// Reads a CSV-file with time series data and fills the container.
    int readDataFromFile(const std::string &file_name);

    std::size_t start_;
    std::size_t end_;
    std::size_t step_size_;
    TimeStepType time_unit_;
    std::vector<std::string> data_unit_string_;
    std::vector<std::size_t> time_steps_;
    std::vector<SensorDataType> vec_names_;
    std::vector< std::vector<float>* > data_vecs_;

};
