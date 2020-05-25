/**
 * \file
 * \author Karsten Rink
 * \date   2013-03-18
 * \brief  Implementation of the StationBorehole class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "BaseLib/Logging.h"

#include "BaseLib/StringTools.h"
#include "BaseLib/DateTools.h"

namespace GeoLib
{

////////////////////////
// The Borehole class //
////////////////////////

StationBorehole::StationBorehole(double x,
                                 double y,
                                 double z,
                                 const std::string& name)
    : Station(x, y, z, name)
{
    type_ = Station::StationType::BOREHOLE;

    // add first point of borehole
    profilePntVec_.push_back(this);
    soilName_.emplace_back("");
}

StationBorehole::~StationBorehole()
{
    // deletes profile vector of borehole, starting at layer 1
    // the first point is NOT deleted as it points to the station object itself
    for (std::size_t k(1); k < profilePntVec_.size(); k++)
    {
        delete profilePntVec_[k];
    }
}

StationBorehole* StationBorehole::createStation(const std::string &line)
{
    StationBorehole* borehole = new StationBorehole();
    std::list<std::string> fields = BaseLib::splitString(line, '\t');

    if (fields.size() >= 5)
    {
        borehole->name_ = fields.front();
        fields.pop_front();
        (*borehole)[0] = strtod(BaseLib::replaceString(",", ".", fields.front()).c_str(), nullptr);
        fields.pop_front();
        (*borehole)[1] = strtod(BaseLib::replaceString(",", ".", fields.front()).c_str(), nullptr);
        fields.pop_front();
        (*borehole)[2] = strtod(BaseLib::replaceString(",", ".", fields.front()).c_str(), nullptr);
        fields.pop_front();
        borehole->depth_ = strtod(BaseLib::replaceString(",", ".", fields.front()).c_str(), nullptr);
        fields.pop_front();
        if (fields.empty())
        {
            borehole->date_ = 0;
        }
        else
        {
            borehole->date_ = BaseLib::strDate2int(fields.front());
            fields.pop_front();
        }
    }
    else
    {
        WARN("Station::createStation() - Unexpected file format.");
        delete borehole;
        return nullptr;
    }
    return borehole;
}

StationBorehole* StationBorehole::createStation(const std::string &name,
                                                double x,
                                                double y,
                                                double z,
                                                double depth,
                                                const std::string &date)
{
    StationBorehole* station = new StationBorehole();
    station->name_  = name;
    (*station)[0]   = x;
    (*station)[1]   = y;
    (*station)[2]   = z;
    station->depth_ = depth;
    if (date != "0000-00-00")
    {
        station->date_ = BaseLib::xmlDate2int(date);
    }
    return station;
}

void StationBorehole::addSoilLayer ( double thickness, const std::string &soil_name)
{
    /*
       // TF - Altmark
       if (profilePntVec_.empty())
        addSoilLayer ((*this)[0], (*this)[1], (*this)[2]-thickness, soil_name);
       else {
        std::size_t idx (profilePntVec_.size());
        // read coordinates from last above
        double x((*profilePntVec_[idx-1])[0]);
        double y((*profilePntVec_[idx-1])[1]);
        double z((*profilePntVec_[idx-1])[2]-thickness);
        addSoilLayer (x, y, z, soil_name);
       }
     */

    // KR - Bode
    if (profilePntVec_.empty())
    {
        addSoilLayer((*this)[0], (*this)[1], (*this)[2], "");
    }

    std::size_t idx (profilePntVec_.size());
    double x((*profilePntVec_[idx - 1])[0]);
    double y((*profilePntVec_[idx - 1])[1]);
    double z((*profilePntVec_[0])[2] - thickness);
    addSoilLayer (x, y, z, soil_name);
}

void StationBorehole::addSoilLayer ( double x, double y, double z, const std::string &soil_name)
{
    profilePntVec_.push_back (new Point (x, y, z));
    soilName_.push_back(soil_name);
}

bool isBorehole(GeoLib::Point const* pnt)
{
    auto const* bh(dynamic_cast<GeoLib::StationBorehole const*>(pnt));
    return bh != nullptr;
}

}  // namespace GeoLib
