/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * \file XmlCndInterface.h
 *
 * Created on 2011-11-23 by Karsten Rink
 */

#ifndef XMLCNDINTERFACE_H
#define XMLCNDINTERFACE_H

#include "XMLInterface.h"

class FEMCondition;

namespace FileIO
{

/**
 * \brief Reads and writes FEM Conditions to and from XML files.
 */
class XmlCndInterface : public XMLInterface
{
public:
	/**
	 * Constructor
	 * \param project Project data.
	 * \param schemaFile An XML schema file (*.xsd) that defines the structure of a valid data file.
	 */
	XmlCndInterface(ProjectData* project, const std::string &schemaFile);

	~XmlCndInterface() {};

	/// Dummy function so class hierarchy works. This needs to be implemented later.
	int readFile(const QString &fileName)
	{
		Q_UNUSED(fileName)
		std::cout << "There is currently no implementation for XmlCndInterface::readFile(const QString&)." << std::endl;
		return 0;
	}

	/// Reads an xml-file containing FEM Conditions such as Boundary- or Initial Conditions
	int readFile(std::vector<FEMCondition*> &conditions, const QString &fileName);

	void setConditionType(FEMCondition::CondType type) { _type = type; };


protected:
	int write(std::ostream& stream);

private:
	/// Read the details of various FEM Conditions from an xml-file
	void readConditions( const QDomNode &condRoot,
	                     std::vector<FEMCondition*> &conditions,
	                     FEMCondition::CondType type);

	QDomElement getCondListElement( QDomDocument doc, QDomElement &root, const QString &condText ) const;
	void writeCondition( QDomDocument doc, QDomElement &listTag, const FEMCondition* cond, const QString &condText, const QString &geoName ) const;

	FEMCondition::CondType _type;
};

}

#endif // XMLCNDINTERFACE_H
