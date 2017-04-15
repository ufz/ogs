/**
 * \file
 * \author Karsten Rink
 * \date   2012-08-01
 * \brief  Implementation of the SensorData class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SensorData.h"

#include <cstdlib>
#include <fstream>

#include <logog/include/logog.hpp>

#include "BaseLib/StringTools.h"
#include "BaseLib/DateTools.h"

SensorData::SensorData(const std::string &file_name)
: _start(0), _end(0), _step_size(0), _time_unit(TimeStepType::NONE)
{
    this->readDataFromFile(file_name);
}

SensorData::SensorData(std::vector<std::size_t> time_steps)
: _start(time_steps[0]), _end(time_steps[time_steps.size()-1]), _step_size(0), _time_unit(TimeStepType::NONE), _time_steps(time_steps)
{
    for (std::size_t i=1; i<time_steps.size(); i++)
    {
        if (time_steps[i-1]>=time_steps[i])
            ERR("Error in SensorData() - Time series has no order!");
    }
}

SensorData::SensorData(std::size_t first_timestep, std::size_t last_timestep, std::size_t step_size)
: _start(first_timestep), _end(last_timestep), _step_size(step_size), _time_unit(TimeStepType::NONE)
{
}

SensorData::~SensorData()
{
    for (std::vector<float>* vec : _data_vecs)
        delete vec;
}


void SensorData::addTimeSeries( const std::string &data_name, std::vector<float> *data, const std::string &data_unit_string )
{
    this->addTimeSeries(SensorData::convertString2SensorDataType(data_name), data, data_unit_string);
}

void SensorData::addTimeSeries(SensorDataType data_name, std::vector<float> *data, const std::string &data_unit_string)
{
    if (_step_size>0) {
        if (((_end-_start)/_step_size) != data->size()) {
            WARN("Warning in SensorData::addTimeSeries() - Lengths of time series does not match number of time steps.");
            return;
        }
    } else {
        if  (data->size() != _time_steps.size()) {
            WARN("Warning in SensorData::addTimeSeries() - Lengths of time series does not match number of time steps.");
            return;
        }
    }

    _vec_names.push_back(data_name);
    _data_vecs.push_back(data);
    _data_unit_string.push_back(data_unit_string);
}

const std::vector<float>* SensorData::getTimeSeries(SensorDataType time_series_name) const
{
    for (std::size_t i=0; i<_vec_names.size(); i++)
    {
        if (time_series_name == _vec_names[i])
            return _data_vecs[i];
    }
    ERR("Error in SensorData::getTimeSeries() - Time series \"%d\" not found.", time_series_name);
    return nullptr;
}

std::string SensorData::getDataUnit(SensorDataType time_series_name) const
{
    for (std::size_t i=0; i<_vec_names.size(); i++)
    {
        if (time_series_name == _vec_names[i])
            return _data_unit_string[i];
    }
    ERR("Error in SensorData::getDataUnit() - Time series \"%d\" not found.", time_series_name);
    return "";
}

int SensorData::readDataFromFile(const std::string &file_name)
{
    std::ifstream in( file_name.c_str() );

    if (!in.is_open())
    {
        INFO("SensorData::readDataFromFile() - Could not open file %s.", file_name.c_str());
        return 0;
    }

    std::string line("");

    /* first line contains field names */
    getline(in, line);
    std::list<std::string> fields = BaseLib::splitString(line, '\t');
    std::list<std::string>::const_iterator it (fields.begin());
    std::size_t nFields = fields.size();

    if (nFields<2)
        return 0;

    std::size_t nDataArrays(nFields-1);

    //create vectors necessary to hold the data
    for (std::size_t i=0; i<nDataArrays; i++)
    {
        this->_vec_names.push_back(SensorData::convertString2SensorDataType(*++it));
        this->_data_unit_string.push_back("");
        auto* data = new std::vector<float>;
        this->_data_vecs.push_back(data);
    }

    while ( getline(in, line) )
    {
        fields = BaseLib::splitString(line, '\t');

        if (nFields == fields.size())
        {
            it = fields.begin();
            std::size_t pos(it->rfind("."));
            std::size_t current_time_step = (pos == std::string::npos) ? atoi((it++)->c_str()) : BaseLib::strDate2int(*it++);
            this->_time_steps.push_back(current_time_step);

            for (std::size_t i=0; i<nDataArrays; i++)
                this->_data_vecs[i]->push_back(
                    static_cast<float>(strtod((it++)->c_str(), nullptr)));
        }
        else
            return 0;
    }

    in.close();

    this->_start = this->_time_steps[0];
    this->_end   = this->_time_steps[this->_time_steps.size()-1];

    return 1;
}

std::string SensorData::convertSensorDataType2String(SensorDataType t)
{
    if (SensorDataType::EVAPORATION == t) return "Evaporation";
    else if (SensorDataType::PRECIPITATION == t) return "Precipitation";
    else if (SensorDataType::TEMPERATURE == t) return "Temperature";
    // pls leave this as last choice
    else return "Unknown";
}

SensorDataType SensorData::convertString2SensorDataType(const std::string &s)
{
    if ((s.compare("Evaporation")==0) || (s.compare("EVAPORATION")==0)) return SensorDataType::EVAPORATION;
    else if ((s.compare("Precipitation")==0) || (s.compare("PRECIPITATION")==0)) return SensorDataType::PRECIPITATION;
    else if ((s.compare("Temperature")==0) || (s.compare("TEMPERATURE")==0)) return SensorDataType::TEMPERATURE;
    else return SensorDataType::OTHER;
}

