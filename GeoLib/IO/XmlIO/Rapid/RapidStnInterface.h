/**
 * \file
 * \author Karsten Rink
 * \date   2012-08-16
 * \brief  Definition of RapidStnInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <vector>

#include "RapidXML/rapidxml.hpp"


namespace GeoLib {
    class Point;
    class StationBorehole;
}

namespace GeoLib
{
namespace IO
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

} // IO
} // GeoLib
