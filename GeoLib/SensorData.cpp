/**
 * \file
 * \author Karsten Rink
 * \date   2012-08-01
 * \brief  Implementation of the SensorData class.
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SensorData.h"

#include <cstdlib>
#include <fstream>

#include "BaseLib/DateTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"

SensorData::SensorData(const std::string& file_name)
    : _start(0), _end(0), _step_size(0), _time_unit(TimeStepType::NONE)
{
    readDataFromFile(file_name);
}

SensorData::SensorData(std::vector<std::size_t> time_steps)
    : _start(time_steps.front()),
      _end(time_steps.back()),
      _step_size(0),
      _time_unit(TimeStepType::NONE),
      _time_steps(time_steps)
{
    if (!std::is_sorted(
            time_steps.begin(), time_steps.end(), std::less_equal{}))
    {
        ERR("Error in SensorData() - Time series has no order!");
    }
}

SensorData::~SensorData()
{
    for (std::vector<float>* vec : _data_vecs)
    {
        delete vec;
    }
}

void SensorData::addTimeSeries(const std::string& data_name,
                               std::vector<float>* data,
                               const std::string& data_unit_string)
{
    addTimeSeries(SensorData::convertString2SensorDataType(data_name),
                  data,
                  data_unit_string);
}

void SensorData::addTimeSeries(SensorDataType data_name,
                               std::vector<float>* data,
                               const std::string& data_unit_string)
{
    if (_step_size > 0)
    {
        if (((_end - _start) / _step_size) != data->size())
        {
            WARN(
                "Warning in SensorData::addTimeSeries() - Lengths of time "
                "series does not match number of time steps.");
            return;
        }
    }
    else
    {
        if (data->size() != _time_steps.size())
        {
            WARN(
                "Warning in SensorData::addTimeSeries() - Lengths of time "
                "series does not match number of time steps.");
            return;
        }
    }

    _vec_names.push_back(data_name);
    _data_vecs.push_back(data);
    _data_unit_string.push_back(data_unit_string);
}

const std::vector<float>* SensorData::getTimeSeries(
    SensorDataType time_series_name) const
{
    for (std::size_t i = 0; i < _vec_names.size(); i++)
    {
        if (time_series_name == _vec_names[i])
        {
            return _data_vecs[i];
        }
    }
    ERR("Error in SensorData::getTimeSeries() - Time series '{:s}' not found.",
        convertSensorDataType2String(time_series_name));
    return nullptr;
}

int SensorData::readDataFromFile(const std::string& file_name)
{
    std::ifstream in(file_name.c_str());

    if (!in.is_open())
    {
        INFO("SensorData::readDataFromFile() - Could not open file {:s}.",
             file_name);
        return 0;
    }

    std::string line;

    /* first line contains field names */
    std::getline(in, line);
    std::list<std::string> fields = BaseLib::splitString(line, '\t');
    std::list<std::string>::const_iterator it(fields.begin());
    std::size_t const nFields = fields.size();

    if (nFields < 2)
    {
        return 0;
    }

    std::size_t const nDataArrays(nFields - 1);

    // create vectors necessary to hold the data
    for (std::size_t i = 0; i < nDataArrays; i++)
    {
        _vec_names.push_back(SensorData::convertString2SensorDataType(*++it));
        _data_unit_string.emplace_back("");
        _data_vecs.push_back(new std::vector<float>);
    }

    while (std::getline(in, line))
    {
        fields = BaseLib::splitString(line, '\t');

        if (nFields != fields.size())
        {
            return 0;
        }

        it = fields.begin();
        std::size_t const pos(it->rfind("."));
        std::size_t const current_time_step = (pos == std::string::npos)
                                                  ? atoi((it++)->c_str())
                                                  : BaseLib::strDate2int(*it++);
        _time_steps.push_back(current_time_step);

        for (std::size_t i = 0; i < nDataArrays; i++)
        {
            _data_vecs[i]->push_back(
                static_cast<float>(strtod((it++)->c_str(), nullptr)));
        }
    }

    in.close();

    _start = _time_steps[0];
    _end = _time_steps[_time_steps.size() - 1];

    return 1;
}

std::string SensorData::convertSensorDataType2String(SensorDataType t)
{
    if (SensorDataType::EVAPORATION == t)
    {
        return "Evaporation";
    }
    if (SensorDataType::PRECIPITATION == t)
    {
        return "Precipitation";
    }
    if (SensorDataType::TEMPERATURE == t)
    {
        return "Temperature";
    }
    // pls leave this as last choice
    return "Unknown";
}

SensorDataType SensorData::convertString2SensorDataType(const std::string& s)
{
    if (s == "Evaporation" || s == "EVAPORATION")
    {
        return SensorDataType::EVAPORATION;
    }
    if (s == "Precipitation" || s == "PRECIPITATION")
    {
        return SensorDataType::PRECIPITATION;
    }
    if (s == "Temperature" || s == "TEMPERATURE")
    {
        return SensorDataType::TEMPERATURE;
    }
    return SensorDataType::OTHER;
}
