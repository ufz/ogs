/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef XMLPRJINTERFACE_H
#define XMLPRJINTERFACE_H

#include <string>

#include <QString>

#include "../XMLInterface.h"
#include "XMLQtInterface.h"

class ProjectData;

namespace FileIO
{

/**
 * \brief Reads and writes project information to and from XML files.
 */
class XmlPrjInterface : public XMLInterface, public XMLQtInterface
{
public:
	/**
	 * Constructor
	 * \param project Project data.
	 */
	XmlPrjInterface(ProjectData &project);

	virtual ~XmlPrjInterface() {}

	/// Reads an xml-file containing a GeoSys project.
	/// Project files currently cover only geo-, msh- and station-data. This will be expanded in the future.
	int readFile(const QString &fileName);

	bool readFile(std::string const& fname) { return readFile(QString(fname.c_str())) != 0; }

	int writeToFile(const std::string& filename);

protected:
	bool write();

private:
	int readInputFiles(QDomNode const& node, QString const& path);

	std::string _filename;

	ProjectData& _project;
};

}

#endif // XMLPRJINTERFACE_H
