/**
 * \file
 * \date Oct 15, 2013
 * \brief Base part of implementation of reading XML files using Qt stuff.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "XMLQtInterface.h"

#include <QByteArray>
#include <QCoreApplication>
#include <QCryptographicHash>
#include <QDir>
#include <QFileInfo>
#include <QXmlSchema>
#include <QXmlSchemaValidator>
#include <QXmlStreamReader>
#include <fstream>
#include <utility>

#include "BaseLib/Logging.h"

namespace BaseLib
{
namespace IO
{
XMLQtInterface::XMLQtInterface(QString schemaFile)
    : schemaFile_(std::move(schemaFile))
{
}

int XMLQtInterface::readFile(const QString& fileName)
{
    fileName_ = fileName;
    QFile file(fileName);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        ERR("XMLQtInterface::readFile(): Can't open xml-file {:s}.",
            fileName.toStdString());
        return 0;
    }
    fileData_ = file.readAll();
    file.close();

    if (!checkHash())
    {
        return 0;
    }

    return 1;
}

int XMLQtInterface::isValid() const
{
    QXmlSchema schema;
    if (schemaFile_.length() > 0)
    {
        auto path =
            QDir(QCoreApplication::applicationDirPath()).filePath(schemaFile_);
        auto url = QUrl::fromLocalFile(path);
        schema.load(url);
    }
    if (schema.isValid())
    {
        QXmlSchemaValidator validator(schema);
        if (validator.validate(fileData_))
        {
            return 1;
        }

        INFO(
            "XMLQtInterface::isValid(): XML file {:s} is invalid (in reference "
            "to schema {:s}).",
            fileName_.toStdString(), schemaFile_.toStdString());
    }
    else
    {
        // The following validator (without constructor arguments) automatically
        // searches for the xsd given in the xml file.
        QXmlSchemaValidator validator;
        if (validator.validate(fileData_))
        {
            return 1;
        }

        INFO(
            "XMLQtInterface::isValid(): XML file {:s} is invalid (in reference "
            "to its schema).",
            fileName_.toStdString());
    }
    return 0;
}

bool XMLQtInterface::checkHash() const
{
    QString md5FileName(fileName_ + ".md5");
    QByteArray fileHash =
        QCryptographicHash::hash(fileData_, QCryptographicHash::Md5);

    QFile file(md5FileName);
    if (file.open(QIODevice::ReadOnly))
    {
        QByteArray referenceHash = file.readAll();
        file.close();
        if (referenceHash == fileHash)
        {
            return true;
        }
        INFO("Hashfile does not match data ... checking file ...");
    }

    if (!this->isValid())
    {
        return false;
    }

    QFile fileMD5(md5FileName);
    if (fileMD5.open(QIODevice::WriteOnly))
    {
        fileMD5.write(fileHash);
        fileMD5.close();
        INFO("File is valid, hashfile written.");
    }
    else
    {
        WARN("File is valid but could not write hashfile!");
    }
    return true;
}

QByteArray const& XMLQtInterface::getContent() const
{
    return fileData_;
}
}  // namespace IO
}  // namespace BaseLib
