/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.com)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.com/LICENSE.txt
 *
 *
 * \file XMLInterface.h
 *
 * Created on 2010-02-18 by Karsten Rink
 */

#ifndef XMLINTERFACE_H
#define XMLINTERFACE_H

#include "ProjectData.h"

#include <QXmlStreamReader>
#include "Writer.h"

class FEMCondition;

class QFile;
class QDomDocument;
class QDomNode;
class QDomElement;

namespace FileIO
{

/**
 * \brief Base class for writing any information to and from XML files.
 */
class XMLInterface : public Writer
{
public:
	/**
	 * Constructor
	 * \param project Project data.
	 * \param schemaFile An XML schema file (*.xsd) that defines the structure of a valid data file.
	 */
	XMLInterface(ProjectData* project, const std::string &schemaFile);

	virtual ~XMLInterface() {};

	/// As QXMLStreamWriter seems currently unable to include style-file links into xml-files, this method will workaround this issue and include the stylefile link.
	int insertStyleFileDefinition(const QString &fileName) const;

	/// Check if the given xml-file is valid considering the schema-file used in the constructor
	int isValid(const QString &fileName) const;

	void setNameForExport(std::string const& name) { _exportName = name; };

	/// Sets the schema filename used to check if xml files are valid.
	void setSchema(const std::string &schemaName);

	/// Reads an xml-file.
	virtual int readFile(const QString &fileName) = 0;

protected:
	/// Checks if a hash for the given data file exists to skip the time-consuming validation part.
	/// If a hash file exists _and_ the hash of the data file is the same as the content of the hash file the validation is skipped
	/// If no hash file exists, the xml-file is validated and a hash file is written if the xml-file was valid.
	bool checkHash(const QString &fileName) const;

	/// Calculates an MD5 hash of the given file.
	QByteArray calcHash(const QString &fileName) const;

	/// Checks if the given file is conform to the given hash.
	bool hashIsGood(const QString &fileName, const QByteArray &hash) const;

	ProjectData* _project;

	std::string _exportName;
	std::string _schemaName;
	std::map<size_t, size_t> _idx_map;
};

}

#endif // XMLINTERFACE_H
