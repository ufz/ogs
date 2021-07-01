/**
 * \file
 * \author Karsten Rink
 * \date   2010-07-01
 * \brief  Implementation of the Station class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
    : Point(x, y, z), _name(std::move(name))
{
}

Station::Station(Point* coords, std::string name)
    : Point(*coords), _name(std::move(name))
{
}

Station::Station(Station const& src)
    : Point(src),
      _name(src._name),
      _station_value(src._station_value),
      _sensor_data(src._sensor_data.get() != nullptr
                       ? new SensorData(*(src._sensor_data.get()))
                       : nullptr)
{
}

Station* Station::createStation(const std::string& line)
{
    std::list<std::string> fields = BaseLib::splitString(line, '\t');

    if (fields.size() < 3)
    {
        INFO("Station::createStation() - Unexpected file format.");
        return nullptr;
    }

    auto it = fields.begin();
    std::string name = *it;
    auto const x = std::strtod(
        (BaseLib::replaceString(",", ".", *(++it))).c_str(), nullptr);
    auto const y = std::strtod(
        (BaseLib::replaceString(",", ".", *(++it))).c_str(), nullptr);
    auto z = 0.0;
    if (++it != fields.end())
    {
        z = std::strtod((BaseLib::replaceString(",", ".", *it)).c_str(),
                        nullptr);
    }
    return new Station(x, y, z, name);
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
