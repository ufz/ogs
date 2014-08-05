/**
 * \file   XmlNumInterface.h
 * \author Karsten Rink
 * \date   2014-08-05
 * \brief  Definition of the XmlNumInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef XMLNUMINTERFACE_H
#define XMLNUMINTERFACE_H

#include "../XMLInterface.h"
#include "XMLQtInterface.h"

namespace FileIO
{

/**
 * \brief Reads and writes GeoObjects to and from XML files.
 */
class XmlNumInterface : public XMLInterface, public XMLQtInterface
{
public:
	/**
	 * Constructor
	 * \param project Project data.
	 */
	XmlNumInterface(ProjectData &project);

	virtual ~XmlNumInterface() {}

	/// Reads an xml-file containing geometric object definitions into the GEOObjects used in the contructor
	int readFile(QString const& fileName);

	bool readFile(std::string const& fname) { return readFile(QString(fname.c_str())) != 0; }

protected:
	void readLinearSolverConfiguration(QDomElement const& lin_root);
	void readIterationScheme(QDomElement const& iteration_root);
	void readConvergenceCriteria(QDomElement const& convergence_root);
	bool write();

private:
	ProjectData &_project;
};

}

#endif // XMLNUMINTERFACE_H
