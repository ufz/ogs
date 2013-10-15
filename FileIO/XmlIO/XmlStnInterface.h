/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Definition of the XmlStnInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef XMLSTNINTERFACE_H
#define XMLSTNINTERFACE_H

#include "RapidXML/rapidxml.hpp"

#include "XMLInterface.h"
#include "XMLQtInterface.h"

namespace GeoLib {
	class StationBorehole;
}

namespace FileIO
{

/**
 * \brief Reads and writes Observation Sites to and from XML files.
 */
class XmlStnInterface : public XMLInterface, public XMLQtInterface
{
public:
	/**
	 * Constructor
	 * \param project Project data.
	 * \param schemaFile An XML schema file (*.xsd) that defines the structure of a valid data file.
	 */
	XmlStnInterface(GeoLib::GEOObjects& geo_objs, const std::string &schemaFile);

	/// Reads an xml-file containing station object definitions into the GEOObjects used in the contructor (requires Qt)
	int readFile(const QString &fileName);

	/// Reads an xml-file using the RapidXML parser integrated in the source code (i.e. this function is usable without Qt)
	int rapidReadFile(const std::string &fileName);

protected:
	int write(std::ostream& stream);

private:
	/// Reads GeoLib::Station- or StationBorehole-objects from an xml-file
	void readStations  ( const QDomNode &stationsRoot, std::vector<GeoLib::Point*>* stations, const std::string &station_file_name);

	/// Writes borehole-specific data to a station-xml-file.
	void writeBoreholeData(QDomDocument &doc,
	                       QDomElement &boreholeTag,
	                       GeoLib::StationBorehole* borehole) const;

	/// Reads the stratigraphy of a borehole from an xml-file
	void readStratigraphy( const QDomNode &stratRoot, GeoLib::StationBorehole*  borehole );

	/// Reads GeoLib::Station- or StationBorehole-objects from an xml-file using the RapidXML parser
	void rapidReadStations(const rapidxml::xml_node<>* station_root, std::vector<GeoLib::Point*> *stations, const std::string &file_name);

	/// Reads the stratigraphy of a borehole from an xml-file using the RapidXML parser
	void rapidReadStratigraphy(const rapidxml::xml_node<>* strat_root, GeoLib::StationBorehole* borehole);

	GeoLib::GEOObjects& _geo_objs;
};

}

#endif // XMLSTNINTERFACE_H
