/**
 * \file
 * \date Oct 15, 2013
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
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
    explicit XMLQtInterface(QString schemaFile = "");
    virtual ~XMLQtInterface() = default;

    /// Reads the file. In an overridden function in the child class be sure to call
    /// XMLQtInterface::readFile(fileName).
    virtual int readFile(const QString &fileName);

protected:
    /// Check if the given xml-file is valid considering the schema-file used in the constructor
    int isValid() const;

    /// Checks if a hash for the given data file exists to skip the time-consuming validation part.
    /// If a hash file exists _and_ the hash of the data file is the same as the content of the hash file the validation is skipped
    /// If no hash file exists, the xml-file is validated and a hash file is written if the xml-file was valid.
    bool checkHash() const;

    QString _schemaFile;

    /// Caches the actual file contents when reading.
    QByteArray _fileData;

    /// The actual file name when reading.
    QString _fileName;
};

} // end namespace IO
} // end namespace BaseLib
