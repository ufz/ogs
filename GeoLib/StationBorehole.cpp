/**
 * \file
 * \author Karsten Rink
 * \date   2013-03-18
 * \brief  Implementation of the StationBorehole class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "StationBorehole.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>

#include "BaseLib/DateTools.h"
#include "BaseLib/Logging.h"
#include "BaseLib/StringTools.h"

namespace GeoLib
{
////////////////////////
// The Borehole class //
////////////////////////

StationBorehole::StationBorehole(double x, double y, double z,
                                 double const depth, const std::string& name,
                                 int date)
    : Station(x, y, z, name), _depth(depth), _date(date)
{
    // add first point of borehole
    _profilePntVec.push_back(this);
    _soilName.emplace_back("");
}

StationBorehole::~StationBorehole()
{
    // deletes profile vector of borehole, starting at layer 1
    // the first point is NOT deleted as it points to the station object itself
    for (std::size_t k(1); k < _profilePntVec.size(); k++)
    {
        delete _profilePntVec[k];
    }
}

StationBorehole* StationBorehole::createStation(const std::string& name,
                                                double x,
                                                double y,
                                                double z,
                                                double depth,
                                                const std::string& date)
{
    int integer_date = 0;
    if (date != "0000-00-00")
    {
        integer_date = BaseLib::xmlDate2int(date);
    }
    return new StationBorehole(x, y, z, depth, name, integer_date);
}

void StationBorehole::addSoilLayer(double thickness,
                                   const std::string& soil_name)
{
    /*
       // TF - Altmark
       if (_profilePntVec.empty())
        addSoilLayer ((*this)[0], (*this)[1], (*this)[2]-thickness, soil_name);
       else {
        std::size_t idx (_profilePntVec.size());
        // read coordinates from last above
        double x((*_profilePntVec[idx-1])[0]);
        double y((*_profilePntVec[idx-1])[1]);
        double z((*_profilePntVec[idx-1])[2]-thickness);
        addSoilLayer (x, y, z, soil_name);
       }
     */

    // KR - Bode
    if (_profilePntVec.empty())
    {
        addSoilLayer((*this)[0], (*this)[1], (*this)[2], soil_name);
    }

    std::size_t idx(_profilePntVec.size());
    double x((*_profilePntVec[idx - 1])[0]);
    double y((*_profilePntVec[idx - 1])[1]);
    double z((*_profilePntVec[0])[2] - thickness);
    addSoilLayer(x, y, z, soil_name);
}

void StationBorehole::addSoilLayer(double x,
                                   double y,
                                   double z,
                                   const std::string& soil_name)
{
    _profilePntVec.push_back(new Point(x, y, z));
    _soilName.push_back(soil_name);
}

bool isBorehole(GeoLib::Point const* pnt)
{
    auto const* bh(dynamic_cast<GeoLib::StationBorehole const*>(pnt));
    return bh != nullptr;
}

}  // namespace GeoLib
