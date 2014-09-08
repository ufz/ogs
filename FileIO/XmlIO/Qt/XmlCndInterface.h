/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Definition of the XmlCndInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef XMLCNDINTERFACE_H
#define XMLCNDINTERFACE_H

// ThirdParty/logog
#include "logog/include/logog.hpp"

// FileIO/XmlIO
#include "../XMLInterface.h"
#include "XMLQtInterface.h"

class FEMCondition;
class ProjectData;

namespace FileIO
{
/**
 * \brief Reads and writes FEM Conditions to and from XML files using the Qt XML parser and validator.
 */
class XmlCndInterface : public XMLInterface, public XMLQtInterface
{
public:
	/**
	 * Constructor
	 * \param project Project data.
	 */
	XmlCndInterface(ProjectData &project);

	~XmlCndInterface() {}

	/// Reads an xml-file containing FEM Conditions such as Boundary- or Initial Conditions
	int readFile(const QString &fileName);

	/// Reads an xml-file containing FEM Conditions (convenience function for using std::strings)
	bool readFile(std::string const& fname) { return readFile(QString(fname.c_str())) != 0; }

	/// Sets the type of condition to be written to a file.
	void setConditionType(FEMCondition::CondType type) { _type = type; }

protected:
	bool write();

private:
	/// Read the details of various FEM Conditions from an xml-file
	void readConditions(const QDomNode &condRoot, FEMCondition::CondType type);

	/// Returns the root node for the kind of FEM Condition specified by condText
	QDomElement getCondListElement(QDomDocument doc, QDomElement &root,
	                               const QString &condText) const;
	
	/// Writes a FEM condition to an xml-file using the write() method
	void writeCondition(QDomDocument doc, QDomElement &listTag, const FEMCondition* cond,
	                    const QString &condText, const QString &geoName) const;

	FEMCondition::CondType _type;

	ProjectData& _project;
};
}

#endif // XMLCNDINTERFACE_H
