/**
 * \file
 * \author Karsten Rink
 * \date   2013-03-18
 * \brief  Implementation of the StationBorehole class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "StationBorehole.h"

#include <cmath>
#include <cstdlib>
#include <fstream>

#include <logog/include/logog.hpp>

#include "BaseLib/StringTools.h"
#include "BaseLib/DateTools.h"

namespace GeoLib
{

////////////////////////
// The Borehole class //
////////////////////////

StationBorehole::StationBorehole(double x, double y, double z, const std::string &name) :
    Station (x, y, z, name), _depth(0), _date(0)
{
    _type = Station::StationType::BOREHOLE;

    // add first point of borehole
    _profilePntVec.push_back(this);
    _soilName.push_back("");
}

StationBorehole::~StationBorehole(void)
{
    // deletes profile vector of borehole, starting at layer 1
    // the first point is NOT deleted as it points to the station object itself
    for (std::size_t k(1); k < _profilePntVec.size(); k++)
        delete _profilePntVec[k];
}

int StationBorehole::find(const std::string &str)
{
    std::size_t size = _soilName.size();
    for (std::size_t i = 0; i < size; i++)
        if (_soilName[i].find(str) == 0)
            return 1;
    return 0;
}

int StationBorehole::readStratigraphyFile(const std::string &path,
                                          std::vector<std::list<std::string> > &data)
{
    std::string line;
    std::ifstream in( path.c_str() );

    if (!in.is_open())
    {
        WARN("StationBorehole::readStratigraphyFile() - Could not open file %s.", path.c_str());
        return 0;
    }

    while ( getline(in, line) )
    {
        std::list<std::string> fields = BaseLib::splitString(line, '\t');
        data.push_back(fields);
    }

    in.close();

    return 1;
}

int StationBorehole::addStratigraphy(const std::string &path, StationBorehole* borehole)
{
    std::vector<std::list<std::string> > data;
    if (readStratigraphyFile(path, data))
    {
        std::size_t size = data.size();
        for (std::size_t i = 0; i < size; i++)
            addLayer(data[i], borehole);

        // check if a layer is missing
        size = borehole->_soilName.size();
        INFO("StationBorehole::addStratigraphy ToDo");
        //    for (std::size_t i=0; i<size; i++)
        //    {
        //        if ((borehole->_soilLayerThickness[i] == -1) ||(borehole->_soilName[i].compare("") == 0))
        //        {
        //            borehole->_soilLayerThickness.clear();
        //            borehole->_soilName.clear();
        //
        //            WARN("StationBorehole::addStratigraphy() - Profile incomplete (Borehole %s, Layer %d missing)", borehole->_name.c_str(), i+1);
        //
        //            return 0;
        //        }
        //    }
    }
    else
        borehole->addSoilLayer(borehole->getDepth(), "depth");

    return 1;
}

int StationBorehole::addLayer(std::list<std::string> fields, StationBorehole* borehole)
{
    if (fields.size() >= 4) /* check if there are enough fields to create a borehole object */
    {
        if (fields.front().compare(borehole->_name) == 0) /* check if the name of the borehole matches the name in the data */
        {
            fields.pop_front();

            // int layer = atoi(fields.front().c_str());
            fields.pop_front();

            ERR("StationBorehole::addLayer - assuming correct order");
            double thickness(strtod(BaseLib::replaceString(",", ".", fields.front()).c_str(), 0));
            fields.pop_front();
            borehole->addSoilLayer(thickness, fields.front());
        }
    }
    else
    {
        WARN("StationBorehole::addLayer() - Unexpected file format (Borehole %s).", borehole->_name.c_str());
        return 0;
    }
    return 1;
}

int StationBorehole::addStratigraphy(const std::vector<Point*> &profile, const std::vector<std::string> &soil_names)
{
    if (((profile.size()-1) == soil_names.size()) && (soil_names.size()>0))
    {
        this->_profilePntVec.push_back(profile[0]);
        std::size_t nLayers = soil_names.size();
        for (std::size_t i=0; i<nLayers; i++)
        {
            this->_profilePntVec.push_back(profile[i+1]);
            this->_soilName.push_back(soil_names[i]);
        }
        return 1;
    }

    ERR("Error in StationBorehole::addStratigraphy() - Length of parameter vectors does not match.");
    return 0;
}

int StationBorehole::addStratigraphies(const std::string &path, std::vector<Point*>* boreholes)
{
    std::vector<std::list<std::string> > data;

    if (readStratigraphyFile(path, data))
    {
        std::string name;

        std::size_t it = 0;
        std::size_t nBoreholes = data.size();
        for (std::size_t i = 0; i < nBoreholes; i++)
        {
            std::list<std::string> fields = data[i];

            if (fields.size() >= 4)
            {
                name = static_cast<StationBorehole*>((*boreholes)[it])->_name;
                if ( fields.front().compare(name) != 0 )
                    if (it < boreholes->size() - 1)
                        it++;

                fields.pop_front();
                //the method just assumes that layers are read in correct order
                fields.pop_front();
                double thickness (strtod(BaseLib::replaceString(",", ".",
                                                       fields.front()).c_str(), 0));
                fields.pop_front();
                std::string soil_name (fields.front());
                fields.pop_front();
                static_cast<StationBorehole*>((*boreholes)[it])->addSoilLayer(
                        thickness,
                        soil_name);
            }
            else
                ERR("Error in StationBorehole::addStratigraphies() - Unexpected file format.");
                //return 0;
        }
    }
    else
        createSurrogateStratigraphies(boreholes);

    return 1;
}

StationBorehole* StationBorehole::createStation(const std::string &line)
{
    StationBorehole* borehole = new StationBorehole();
    std::list<std::string> fields = BaseLib::splitString(line, '\t');

    if (fields.size() >= 5)
    {
        borehole->_name = fields.front();
        fields.pop_front();
        (*borehole)[0] = strtod(BaseLib::replaceString(",", ".", fields.front()).c_str(), nullptr);
        fields.pop_front();
        (*borehole)[1] = strtod(BaseLib::replaceString(",", ".", fields.front()).c_str(), nullptr);
        fields.pop_front();
        (*borehole)[2] = strtod(BaseLib::replaceString(",", ".", fields.front()).c_str(), nullptr);
        fields.pop_front();
        borehole->_depth = strtod(BaseLib::replaceString(",", ".", fields.front()).c_str(), nullptr);
        fields.pop_front();
        if (fields.empty())
            borehole->_date = 0;
        else
        {
            borehole->_date = BaseLib::strDate2int(fields.front());
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
    station->_name  = name;
    (*station)[0]   = x;
    (*station)[1]   = y;
    (*station)[2]   = z;
    station->_depth = depth;
    if (date.compare("0000-00-00") != 0)
        station->_date  = BaseLib::xmlDate2int(date);
    return station;
}

void StationBorehole::createSurrogateStratigraphies(std::vector<Point*>* boreholes)
{
    std::size_t nBoreholes = boreholes->size();
    for (std::size_t i = 0; i < nBoreholes; i++)
    {
        StationBorehole* bore = static_cast<StationBorehole*>((*boreholes)[i]);
        bore->addSoilLayer(bore->getDepth(), "depth");
    }
}

void StationBorehole::addSoilLayer ( double thickness, const std::string &soil_name)
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
        addSoilLayer ((*this)[0], (*this)[1], (*this)[2], "");

    std::size_t idx (_profilePntVec.size());
    double x((*_profilePntVec[idx - 1])[0]);
    double y((*_profilePntVec[idx - 1])[1]);
    double z((*_profilePntVec[0])[2] - thickness);
    addSoilLayer (x, y, z, soil_name);
}

void StationBorehole::addSoilLayer ( double x, double y, double z, const std::string &soil_name)
{
    _profilePntVec.push_back (new Point (x, y, z));
    _soilName.push_back(soil_name);
}

bool isBorehole(GeoLib::Point const* pnt)
{
    GeoLib::StationBorehole const* bh(
        dynamic_cast<GeoLib::StationBorehole const*>(pnt));
    return bh != nullptr;
}

} // namespace
