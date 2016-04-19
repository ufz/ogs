/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Definition of the XmlGspInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef XMLGSPINTERFACE_H
#define XMLGSPINTERFACE_H

#include <functional>
#include <vector>
#include <string>

#include <QString>

#include "../XMLInterface.h"
#include "XMLQtInterface.h"

namespace MeshLib
{
class Mesh;
}

namespace GeoLib
{
class GEOObjects;
}

namespace FileIO
{

/**
 * \brief Reads and writes project information to and from XML files.
 */
class XmlGspInterface : public XMLInterface, public XMLQtInterface
{
public:
	XmlGspInterface(
	    GeoLib::GEOObjects& geoObjects,
	    std::vector<MeshLib::Mesh*> const& mesh_vector,
	    std::function<void(MeshLib::Mesh* const)>&& add_mesh_callback);

	virtual ~XmlGspInterface() {}

	/// Reads an xml-file containing a GeoSys project.
	/// Project files currently cover only geo-, msh- and station-data. This will be expanded in the future.
	int readFile(const QString &fileName);

	bool readFile(std::string const& fname) { return readFile(QString(fname.c_str())) != 0; }

	int writeToFile(const std::string& filename);

protected:
	bool write();

private:
	std::string _filename;
	GeoLib::GEOObjects& _geoObjects;
	std::vector<MeshLib::Mesh*> const& _mesh_vector;
	std::function<void(MeshLib::Mesh* const)> _add_mesh_callback;
};

}

#endif // XMLGSPINTERFACE_H
