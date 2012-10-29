/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file RapidStnInterface.h
 *
 *  Created on 2012-08-16 by Karsten Rink
 */

#ifndef RAPIDXMLINTERFACE_H
#define RAPIDXMLINTERFACE_H

#include "RapidXML/rapidxml.hpp"
#include <vector>

#include "Point.h"

namespace GeoLib {
	class StationBorehole;
}

namespace FileIO
{

/**
 * \brief Base class for writing any information to and from XML files.
 */
class RapidStnInterface
{
public:
	/// Reads an xml-file using the RapidXML parser integrated in the source code (i.e. this function is usable without Qt)
	//int rapidReadFile(const std::string &fileName);
	static std::vector<GeoLib::Point*> *readStationFile(const std::string &fileName);

private:
	/// Reads GEOLIB::Station- or StationBorehole-objects from an xml-file using the RapidXML parser
	static void readStations(const rapidxml::xml_node<>* station_root, std::vector<GeoLib::Point*> *stations, const std::string &file_name);

	/// Reads the stratigraphy of a borehole from an xml-file using the RapidXML parser
	static void readStratigraphy(const rapidxml::xml_node<>* strat_root, GeoLib::StationBorehole* borehole);
};

}

#endif // RAPIDXMLINTERFACE_H
