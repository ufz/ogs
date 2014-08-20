/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Definition of the XmlProjectInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef XMLPROJECTINTERFACE_H
#define XMLPROJECTINTERFACE_H

#include "../XMLInterface.h"
#include "XMLQtInterface.h"

namespace FileIO
{

/**
 * \brief Reads and writes project information to and from XML files.
 */
class XmlProjectInterface : public XMLInterface, public XMLQtInterface
{
public:
	/**
	 * Constructor
	 * \param project Project data.
	 */
	XmlProjectInterface(ProjectData &project);

	virtual ~XmlProjectInterface() {};

	/// Reads an xml-file containing a GeoSys project.
	/// Project files currently cover only geo-, msh- and station-data. This will be expanded in the future.
	int readFile(const QString &fileName);

	bool readFile(std::string const& fname) { return readFile(QString(fname.c_str())) != 0; }

	int writeToFile(std::string filename);

protected:
	bool write();

	// helper function for writing file + filename of every individual file-type
	void writeFile(QDomDocument &doc,
	               QDomElement &parent_node,
	               QString const& tag_name,
	               XMLInterface &xml_interface, 
	               std::string const& path, 
	               std::string const& file_name, 
	               std::string const& extension);

private:
	std::string _filename;

	ProjectData& _project;
};

}

#endif // XMLPROJECTINTERFACE_H
