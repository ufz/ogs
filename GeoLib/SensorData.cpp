/**
 * \file
 * \author Karsten Rink
 * \date   2012-08-01
 * \brief  Implementation of the SensorData class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SensorData.h"

#include <cstdlib>
#include <fstream>

#include "BaseLib/Logging.h"

#include "BaseLib/StringTools.h"
#include "BaseLib/DateTools.h"

SensorData::SensorData(const std::string &file_name)
: start_(0), end_(0), step_size_(0), time_unit_(TimeStepType::NONE)
{
    this->readDataFromFile(file_name);
}

SensorData::SensorData(std::vector<std::size_t> time_steps)
: start_(time_steps[0]), end_(time_steps[time_steps.size()-1]), step_size_(0), time_unit_(TimeStepType::NONE), time_steps_(time_steps)
{
    for (std::size_t i=1; i<time_steps.size(); i++)
    {
        if (time_steps[i - 1] >= time_steps[i])
        {
            ERR("Error in SensorData() - Time series has no order!");
        }
    }
}

SensorData::SensorData(std::size_t first_timestep, std::size_t last_timestep, std::size_t step_size)
: start_(first_timestep), end_(last_timestep), step_size_(step_size), time_unit_(TimeStepType::NONE)
{
}

SensorData::~SensorData()
{
    for (std::vector<float>* vec : data_vecs_)
    {
        delete vec;
    }
}


void SensorData::addTimeSeries( const std::string &data_name, std::vector<float> *data, const std::string &data_unit_string )
{
    this->addTimeSeries(SensorData::convertString2SensorDataType(data_name), data, data_unit_string);
}

void SensorData::addTimeSeries(SensorDataType data_name, std::vector<float> *data, const std::string &data_unit_string)
{
    if (step_size_>0) {
        if (((end_-start_)/step_size_) != data->size()) {
            WARN("Warning in SensorData::addTimeSeries() - Lengths of time series does not match number of time steps.");
            return;
        }
    } else {
        if  (data->size() != time_steps_.size()) {
            WARN("Warning in SensorData::addTimeSeries() - Lengths of time series does not match number of time steps.");
            return;
        }
    }

    vec_names_.push_back(data_name);
    data_vecs_.push_back(data);
    data_unit_string_.push_back(data_unit_string);
}

const std::vector<float>* SensorData::getTimeSeries(SensorDataType time_series_name) const
{
    for (std::size_t i=0; i<vec_names_.size(); i++)
    {
        if (time_series_name == vec_names_[i])
        {
            return data_vecs_[i];
        }
    }
    ERR("Error in SensorData::getTimeSeries() - Time series '{:s}' not found.",
        convertSensorDataType2String(time_series_name));
    return nullptr;
}

int SensorData::readDataFromFile(const std::string &file_name)
{
    std::ifstream in( file_name.c_str() );

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
    std::list<std::string>::const_iterator it (fields.begin());
    std::size_t nFields = fields.size();

    if (nFields < 2)
    {
        return 0;
    }

    std::size_t nDataArrays(nFields-1);

    //create vectors necessary to hold the data
    for (std::size_t i=0; i<nDataArrays; i++)
    {
        this->vec_names_.push_back(SensorData::convertString2SensorDataType(*++it));
        this->data_unit_string_.emplace_back("");
        auto* data = new std::vector<float>;
        this->data_vecs_.push_back(data);
    }

    while (std::getline(in, line))
    {
        fields = BaseLib::splitString(line, '\t');

        if (nFields == fields.size())
        {
            it = fields.begin();
            std::size_t pos(it->rfind("."));
            std::size_t current_time_step = (pos == std::string::npos) ? atoi((it++)->c_str()) : BaseLib::strDate2int(*it++);
            this->time_steps_.push_back(current_time_step);

            for (std::size_t i = 0; i < nDataArrays; i++)
            {
                this->data_vecs_[i]->push_back(
                    static_cast<float>(strtod((it++)->c_str(), nullptr)));
            }
        }
        else
        {
            return 0;
        }
    }

    in.close();

    this->start_ = this->time_steps_[0];
    this->end_   = this->time_steps_[this->time_steps_.size()-1];

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

SensorDataType SensorData::convertString2SensorDataType(const std::string &s)
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

