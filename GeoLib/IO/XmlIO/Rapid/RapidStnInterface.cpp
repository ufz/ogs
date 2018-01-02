/**
 * \file
 * \author Karsten Rink
 * \date   2012-08-16
 * \brief  Implementation of the RapidStnInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RapidStnInterface.h"

#include <fstream>

#include "BaseLib/StringTools.h"
#include "BaseLib/FileTools.h"

#include "GeoLib/Station.h"
#include "GeoLib/StationBorehole.h"

namespace GeoLib
{
namespace IO
{

std::vector<GeoLib::Point*> *RapidStnInterface::readStationFile(const std::string &fileName)
{
    auto* stations = new std::vector<GeoLib::Point*>;
    std::ifstream in(fileName.c_str());
    if (in.fail())
    {
        ERR("XmlStnInterface::rapidReadFile() - Can't open xml-file.");
        return nullptr;
    }

    // read the file in a buffer
    std::stringstream sstr;
    sstr << in.rdbuf();
    std::string buffer = sstr.str();
    in.close();

    // build DOM tree
    rapidxml::xml_document<> doc;
    doc.parse<rapidxml::parse_non_destructive>(
        const_cast<char*>(buffer.data()));

    // parse content
    if (std::string(doc.first_node()->name()) != "OpenGeoSysSTN")
    {
        ERR("XmlStnInterface::readFile() - Unexpected XML root.");
        return nullptr;
    }

    // iterate over all station lists
    for (rapidxml::xml_node<>* station_list = doc.first_node()->first_node(); station_list; station_list = station_list->next_sibling())
    {
        std::string stnName("[NN]");

        stnName = station_list->first_node("name")->value();
        for (rapidxml::xml_node<>* list_item = station_list->first_node(); list_item; list_item = list_item->next_sibling())
        {
            std::string b(list_item->name());
            if (b == "stations")
                RapidStnInterface::readStations(list_item, stations, fileName);
            if (b == "boreholes")
                RapidStnInterface::readStations(list_item, stations, fileName);
        }
    }

    doc.clear();

    return stations;
}
/*
int RapidStnInterface::rapidReadFile(const std::string &fileName)
{
    GEOLIB::GEOObjects* geoObjects = _project->getGEOObjects();

    std::ifstream in(fileName.c_str());
    if (in.fail())
    {
        ERR("XmlStnInterface::rapidReadFile() - Can't open xml-file.");
        return 0;
    }

    // buffer file
    in.seekg(0, std::ios::end);
    std::size_t length = in.tellg();
    in.seekg(0, std::ios::beg);
    char* buffer = new char[length+1];
    in.read(buffer, length);
    buffer[in.gcount()] = '\0';
    in.close();

    // build DOM tree
    rapidxml::xml_document<> doc;
    doc.parse<0>(buffer);

    // parse content
    if (std::string(doc.first_node()->name()).compare("OpenGeoSysSTN"))
    {
        ERR("XmlStnInterface::readFile() - Unexpected XML root.");
        return 0;
    }

    // iterate over all station lists
    for (rapidxml::xml_node<>* station_list = doc.first_node()->first_node(); station_list; station_list = station_list->next_sibling())
    {
        std::vector<GEOLIB::Point*>* stations = new std::vector<GEOLIB::Point*>;
        std::string stnName("[NN]");

        stnName = station_list->first_node("name")->value();
        for (rapidxml::xml_node<>* list_item = station_list->first_node(); list_item; list_item = list_item->next_sibling())
        {
            std::string b(list_item->name());
            if (std::string(list_item->name()).compare("stations") == 0)
                XmlStnInterface::rapidReadStations(list_item, stations, fileName);
            if (std::string(list_item->name()).compare("boreholes") == 0)
                XmlStnInterface::rapidReadStations(list_item, stations, fileName);
        }

        if (!stations->empty())
            geoObjects->addStationVec(stations, stnName);
        else
            delete stations;
    }

    doc.clear();
    delete [] buffer;

    return 1;
}
*/
void RapidStnInterface::readStations(const rapidxml::xml_node<>* station_root, std::vector<GeoLib::Point*> *stations, const std::string &file_name)
{
    for (rapidxml::xml_node<>* station_node = station_root->first_node(); station_node; station_node = station_node->next_sibling())
    {
        if (station_node->first_attribute("id") && station_node->first_attribute("x") && station_node->first_attribute("y"))
        {
            double zVal(0.0);
            if (station_node->first_attribute("z"))
                zVal = strtod(station_node->first_attribute("z")->value(),
                              nullptr);

            std::string station_name(""), sensor_data_file_name(""), bdate_str("0000-00-00");
            double station_value(0.0), borehole_depth(0.0);
            if (station_node->first_node("name"))
                station_name = station_node->first_node("name")->value();
            if (station_node->first_node("sensordata"))
                sensor_data_file_name = station_node->first_node("sensordata")->value();
            if (station_node->first_node("value"))
                station_value =
                    strtod(station_node->first_node("value")->value(), nullptr);
            /* add other station features here */

            if (std::string(station_node->name()) == "station")
            {
                auto* s = new GeoLib::Station(
                    strtod(station_node->first_attribute("x")->value(),
                           nullptr),
                    strtod(station_node->first_attribute("y")->value(),
                           nullptr),
                    zVal,
                    station_name);
                s->setStationValue(station_value);
                if (!sensor_data_file_name.empty())
                    s->addSensorDataFromCSV(BaseLib::copyPathToFileName(sensor_data_file_name, file_name));
                stations->push_back(s);
            }
            else if (std::string(station_node->name()) == "borehole")
            {
                if (station_node->first_node("bdepth"))
                    borehole_depth = strtod(
                        station_node->first_node("bdepth")->value(), nullptr);
                if (station_node->first_node("bdate"))
                    bdate_str = station_node->first_node("bdate")->value();
                /* add other borehole features here */

                GeoLib::StationBorehole* s =
                    GeoLib::StationBorehole::createStation(
                        station_name,
                        strtod(station_node->first_attribute("x")->value(),
                               nullptr),
                        strtod(station_node->first_attribute("y")->value(),
                               nullptr),
                        zVal,
                        borehole_depth,
                        bdate_str);
                s->setStationValue(station_value);

                if (station_node->first_node("strat"))
                    RapidStnInterface::readStratigraphy(station_node->first_node("strat"), s);

                stations->push_back(s);

            }
        }
        else
        {
            ERR(
                "XmlStnInterface::rapidReadStations() - Attribute missing in "
                "<station> tag ...");
        }
    }
}

void RapidStnInterface::readStratigraphy( const rapidxml::xml_node<>* strat_root, GeoLib::StationBorehole* borehole )
{
    double depth_check((*borehole)[2]);

    for (rapidxml::xml_node<>* horizon_node = strat_root->first_node("horizon"); horizon_node; horizon_node = horizon_node->next_sibling())
    {
        if (horizon_node->first_attribute("id") && horizon_node->first_attribute("x") &&
            horizon_node->first_attribute("y")  && horizon_node->first_attribute("z"))
        {
            std::string horizon_name("[NN]");
            if (horizon_node->first_node("name"))
                horizon_name = horizon_node->first_node("name")->value();
            /* add other horizon features here */

            double depth(
                strtod(horizon_node->first_attribute("z")->value(), nullptr));
            if (fabs(depth - depth_check) > std::numeric_limits<double>::epsilon()) // skip soil-layer if its thickness is zero
            {
                borehole->addSoilLayer(
                    strtod(horizon_node->first_attribute("x")->value(),
                           nullptr),
                    strtod(horizon_node->first_attribute("y")->value(),
                           nullptr),
                    depth,
                    horizon_name);
                depth_check = depth;
            }
            else
            {
                WARN(
                    "Warning: Skipped layer \"%s\" in borehole \"%s\" because "
                    "of thickness 0.0.",
                    horizon_name.c_str(), borehole->getName().c_str());
            }
        }
        else
        {
            WARN(
                "XmlStnInterface::rapidReadStratigraphy() - Attribute missing "
                "in <horizon> tag ...");
        }
    }
}

} // IO
} // GeoLib
