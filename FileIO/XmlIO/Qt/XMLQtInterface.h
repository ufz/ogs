/**
 * @file
 * @author git blame XMLQtInterface.h
 * @date Oct 15, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */
#ifndef XMLQTINTERFACE_H_
#define XMLQTINTERFACE_H_

#include <string>

#include <QString>
#include <QByteArray>

class QFile;
class QDomDocument;
class QDomNode;
class QDomElement;

namespace FileIO
{

class XMLQtInterface
{
public:
	XMLQtInterface(const std::string &schemaFile = "");
	virtual ~XMLQtInterface() {}

	/// As QXMLStreamWriter seems currently unable to include style-file links into xml-files, this method will workaround this issue and include the stylefile link.
	int insertStyleFileDefinition(const QString &fileName) const;

	/// Check if the given xml-file is valid considering the schema-file used in the constructor
	int isValid(const QString &fileName) const;

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
	bool isHashGood(const QString &fileName, const QByteArray &hash) const;

	std::string _schemaName;
};

} // end namespace FileIO

#endif /* XMLQTINTERFACE_H_ */
