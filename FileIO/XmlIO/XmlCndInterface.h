/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Definition of the XmlCndInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "XMLInterface.h"
#include "XMLQtInterface.h"

class FEMCondition;
class ProjectData;

namespace FileIO
{
/**
 * \brief Reads and writes FEM Conditions to and from XML files.
 */
class XmlCndInterface : public XMLInterface, public XMLQtInterface
{
public:
	/**
	 * Constructor
	 * \param project Project data.
	 * \param schemaFile An XML schema file (*.xsd) that defines the structure of a valid data file.
	 */
	XmlCndInterface(ProjectData* project, const std::string &schemaFile);

	~XmlCndInterface() {}

	/// Dummy function so class hierarchy works. This needs to be implemented later.
	int readFile(const QString &fileName)
	{
		Q_UNUSED(fileName)
		INFO("There is currently no implementation for XmlCndInterface::readFile(const QString&).");
		return 0;
	}

	/// Reads an xml-file containing FEM Conditions such as Boundary- or Initial Conditions
	int readFile(const QString &fileName, std::vector<FEMCondition*> &conditions);

	void setConditionType(FEMCondition::CondType type) { _type = type; };


protected:
	int write(std::ostream& stream);

private:
	/// Read the details of various FEM Conditions from an xml-file
	void readConditions(const QDomNode &condRoot, std::vector<FEMCondition*> &conditions,
	                    FEMCondition::CondType type);

	QDomElement getCondListElement(QDomDocument doc, QDomElement &root,
	                               const QString &condText) const;
	void writeCondition(QDomDocument doc, QDomElement &listTag, const FEMCondition* cond,
	                    const QString &condText, const QString &geoName) const;

	FEMCondition::CondType _type;

	ProjectData* _project;
};
}

#endif // XMLCNDINTERFACE_H
