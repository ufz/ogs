/**
 * \file
 * \author Karsten Rink
 * \date   2010-07-01
 * \brief  Implementation of the Station class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Station.h"

#include <cstdlib>
#include <utility>

#include "BaseLib/Logging.h"

#include "BaseLib/StringTools.h"

namespace GeoLib
{
Station::Station(double x, double y, double z, std::string name)
    : Point(x, y, z), name_(std::move(name))

{}

Station::Station(Point* coords, std::string name)
    : Point(*coords), name_(std::move(name))
{}

Station::Station(Station const& src)
    : Point(src),
      name_(src.name_),
      type_(src.type_),
      station_value_(src.station_value_)
{}

Station::~Station()
{
    delete this->sensor_data_;
}

Station* Station::createStation(const std::string & line)
{
    Station* station = new Station();
    std::list<std::string> fields = BaseLib::splitString(line, '\t');

    if (fields.size() >= 3)
    {
        auto it = fields.begin();
        station->name_ = *it;
        (*station)[0] = std::strtod((BaseLib::replaceString(",", ".", *(++it))).c_str(), nullptr);
        (*station)[1] = std::strtod((BaseLib::replaceString(",", ".", *(++it))).c_str(), nullptr);
        if (++it != fields.end())
        {
            (*station)[2] = std::strtod(
                (BaseLib::replaceString(",", ".", *it)).c_str(), nullptr);
        }
    }
    else
    {
        INFO("Station::createStation() - Unexpected file format.");
        delete station;
        return nullptr;
    }
    return station;
}

Station* Station::createStation(const std::string& name, double x, double y,
                                double z)
{
    return new Station(x, y, z, name);
}

bool isStation(GeoLib::Point const* pnt)
{
    auto const* bh(dynamic_cast<GeoLib::Station const*>(pnt));
    return bh != nullptr;
}

}  // namespace GeoLib
