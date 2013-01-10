/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Definition of the XmlGspInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef XMLGSPINTERFACE_H
#define XMLGSPINTERFACE_H

#include "XMLInterface.h"

namespace FileIO
{

/**
 * \brief Reads and writes project information to and from XML files.
 */
class XmlGspInterface : public XMLInterface
{
public:
	/**
	 * Constructor
	 * \param project Project data.
	 * \param schemaFile An XML schema file (*.xsd) that defines the structure of a valid data file.
	 */
	XmlGspInterface(ProjectData* project, const std::string &schemaFile);

	virtual ~XmlGspInterface() {};

	/// Reads an xml-file containing a GeoSys project.
	/// Project files currently cover only geo-, msh- and station-data. This will be expanded in the future.
	int readFile(const QString &fileName);

	int writeToFile(std::string filename);

protected:
	int write(std::ostream& stream);

private:
	std::string _filename;
};

}

#endif // XMLGSPINTERFACE_H
