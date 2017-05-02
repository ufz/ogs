/**
 * @file
 * @author git blame XMLQtInterface.h
 * @date Oct 15, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */
#pragma once

#include <string>

#include <QString>
#include <QByteArray>

class QFile;
class QDomDocument;
class QDomNode;
class QDomElement;

namespace BaseLib
{
namespace IO
{

class XMLQtInterface
{
public:
    XMLQtInterface(std::string schemaFile = "");
    virtual ~XMLQtInterface() = default;

    /// As QXMLStreamWriter seems currently unable to include style-file links into xml-files, this method will workaround this issue and include the stylefile link.
    int insertStyleFileDefinition(const QString &fileName) const;

    /// Reads the file. In an overriden function in the child class be sure to call
    /// XMLQtInterface::readFile(fileName).
    virtual int readFile(const QString &fileName);

protected:
    /// Check if the given xml-file is valid considering the schema-file used in the constructor
    int isValid() const;

    /// Checks if a hash for the given data file exists to skip the time-consuming validation part.
    /// If a hash file exists _and_ the hash of the data file is the same as the content of the hash file the validation is skipped
    /// If no hash file exists, the xml-file is validated and a hash file is written if the xml-file was valid.
    bool checkHash() const;

    /// Checks if the given file is conform to the given hash.
    bool isHashGood(const QByteArray &hash) const;

    std::string _schemaName;

    /// Caches the actual file contents when reading.
    QByteArray _fileData;

    /// The actual file name when reading.
    QString _fileName;
};

} // end namespace IO
} // end namespace BaseLib
