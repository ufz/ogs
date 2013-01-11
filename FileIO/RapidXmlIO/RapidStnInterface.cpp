/**
 * \file
 * \author Karsten Rink
 * \date   2012-08-16
 * \brief  Implementation of the RapidStnInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <iostream>

//RapidXML
#include "RapidStnInterface.h"

// BaseLib
#include "StringTools.h"
#include "FileTools.h"

// GeoLib
#include "Station.h"

namespace FileIO
{

std::vector<GeoLib::Point*> *RapidStnInterface::readStationFile(const std::string &fileName)
{
	std::vector<GeoLib::Point*> *stations = new std::vector<GeoLib::Point*>;
	std::ifstream in(fileName.c_str());
	if (in.fail())
	{
		std::cout << "XmlStnInterface::rapidReadFile() - Can't open xml-file." << std::endl;
		return NULL;
	}

	// buffer file
	in.seekg(0, std::ios::end);
	size_t length = in.tellg();
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
		std::cout << "XmlStnInterface::readFile() - Unexpected XML root." << std::endl;
		return NULL;
	}

	// iterate over all station lists
	for (rapidxml::xml_node<>* station_list = doc.first_node()->first_node(); station_list; station_list = station_list->next_sibling())
	{
		std::string stnName("[NN]");

		stnName = station_list->first_node("name")->value();
		for (rapidxml::xml_node<>* list_item = station_list->first_node(); list_item; list_item = list_item->next_sibling())
		{
			std::string b(list_item->name());
			if (std::string(list_item->name()).compare("stations") == 0)
				RapidStnInterface::readStations(list_item, stations, fileName);
			if (std::string(list_item->name()).compare("boreholes") == 0)
				RapidStnInterface::readStations(list_item, stations, fileName);
		}
	}

	doc.clear();
	delete [] buffer;

	return stations;
}
/*
int RapidStnInterface::rapidReadFile(const std::string &fileName)
{
	GEOLIB::GEOObjects* geoObjects = _project->getGEOObjects();

	std::ifstream in(fileName.c_str());
	if (in.fail())
	{
		std::cout << "XmlStnInterface::rapidReadFile() - Can't open xml-file." << std::endl;
		return 0;
	}

	// buffer file
	in.seekg(0, std::ios::end);
	size_t length = in.tellg();
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
		std::cout << "XmlStnInterface::readFile() - Unexpected XML root." << std::endl;
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
				zVal = strtod(station_node->first_attribute("z")->value(), 0);

			std::string station_name(""), sensor_data_file_name(""), bdate_str("0000-00-00");
			double station_value(0.0), borehole_depth(0.0);
			if (station_node->first_node("name"))
				station_name = station_node->first_node("name")->value();
			if (station_node->first_node("sensordata"))
				sensor_data_file_name = station_node->first_node("sensordata")->value();
			if (station_node->first_node("value"))
				station_value = strtod(station_node->first_node("value")->value(), 0);
			/* add other station features here */

			if (std::string(station_node->name()).compare("station") == 0)
			{
				GeoLib::Station* s = new GeoLib::Station(
							strtod(station_node->first_attribute("x")->value(), 0),
				            strtod(station_node->first_attribute("y")->value(), 0),
				            zVal,
							station_name);
				s->setStationValue(station_value);
				if (!sensor_data_file_name.empty())
					s->addSensorDataFromCSV(BaseLib::copyPathToFileName(sensor_data_file_name, file_name));
				stations->push_back(s);
			}
			else if (std::string(station_node->name()).compare("borehole") == 0)
			{
				if (station_node->first_node("bdepth"))
					borehole_depth = strtod(station_node->first_node("bdepth")->value(), 0);
				if (station_node->first_node("bdate"))
					bdate_str = station_node->first_node("bdate")->value();
				/* add other borehole features here */

				GeoLib::StationBorehole* s = GeoLib::StationBorehole::createStation(
				        station_name,
				        strtod(station_node->first_attribute("x")->value(), 0),
				        strtod(station_node->first_attribute("y")->value(), 0),
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
			std::cout << "XmlStnInterface::rapidReadStations() - Attribute missing in <station> tag ..." << std::endl;
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

			double depth (strtod(horizon_node->first_attribute("z")->value(), 0));
			if (fabs(depth - depth_check) > std::numeric_limits<double>::epsilon()) // skip soil-layer if its thickness is zero
			{
				borehole->addSoilLayer(strtod(horizon_node->first_attribute("x")->value(), 0),
									   strtod(horizon_node->first_attribute("y")->value(), 0),
									   depth,
									   horizon_name);
				depth_check = depth;
			}
			else
				std::cout << "Warning: Skipped layer \"" << horizon_name << "\" in borehole \""
					      << borehole->getName() << "\" because of thickness 0.0." << std::endl;
		}
		else
			std::cout <<
			"XmlStnInterface::rapidReadStratigraphy() - Attribute missing in <horizon> tag ..." << std::endl;
	}
}

}
