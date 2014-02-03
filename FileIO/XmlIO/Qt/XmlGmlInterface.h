/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Definition of the XmlGmlInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef XMLGMLINTERFACE_H
#define XMLGMLINTERFACE_H

#include "../XMLInterface.h"
#include "XMLQtInterface.h"

namespace FileIO
{

/**
 * \brief Reads and writes GeoObjects to and from XML files.
 */
class XmlGmlInterface : public XMLInterface, public XMLQtInterface
{
public:
	/**
	 * Constructor
	 * \param project Project data.
	 */
	XmlGmlInterface(GeoLib::GEOObjects& geo_objs);

	virtual ~XmlGmlInterface() {}

	/// Reads an xml-file containing geometric object definitions into the GEOObjects used in the contructor
	int readFile(const QString &fileName);

	bool readFile(std::string const& fname) { return readFile(QString(fname.c_str())) != 0; }

protected:
	int write(std::ostream& stream);

private:
	/// Reads GeoLib::Point-objects from an xml-file
	void readPoints    ( const QDomNode &pointsRoot,
	                     std::vector<GeoLib::Point*>* points,
	                     std::map<std::string, std::size_t>* &pnt_names );

	/// Reads GeoLib::Polyline-objects from an xml-file
	void readPolylines ( const QDomNode &polylinesRoot,
	                     std::vector<GeoLib::Polyline*>* polylines,
	                     std::vector<GeoLib::Point*>* points,
	                     const std::vector<std::size_t> &pnt_id_map,
	                     std::map<std::string, std::size_t>* &ply_names );

	/// Reads GeoLib::Surface-objects from an xml-file
	void readSurfaces  ( const QDomNode &surfacesRoot,
	                     std::vector<GeoLib::Surface*>* surfaces,
	                     std::vector<GeoLib::Point*>* points,
	                     const std::vector<std::size_t> &pnt_id_map,
	                     std::map<std::string, std::size_t>* &sfc_names );

	GeoLib::GEOObjects& _geo_objs;
	std::map<std::size_t, std::size_t> _idx_map;
};

}

#endif // XMLGMLINTERFACE_H
