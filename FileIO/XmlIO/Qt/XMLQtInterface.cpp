/**
 * @file
 * @author git blame XMLQtInterface.cpp
 * @date Oct 15, 2013
 * @brief Base part of implementation of reading XML files using Qt stuff.
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "XmlIO/Qt/XMLQtInterface.h"

#include <fstream>

#include <QXmlStreamReader>
#include <QXmlSchema>
#include <QXmlSchemaValidator>
#include <QFileInfo>
#include <QByteArray>
#include <QCryptographicHash>

// ThirdParty/logog
#include "logog/include/logog.hpp"

namespace FileIO
{

XMLQtInterface::XMLQtInterface(const std::string &schemaFile) :
		_schemaName(schemaFile)
{}

int XMLQtInterface::readFile(const QString &fileName)
{
	_fileName = fileName;
	QFile file(fileName);
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
	{
		ERR("XMLQtInterface::readFile(): Can't open xml-file %s.", fileName.data());
		return 0;
	}
	_fileData = file.readAll();

	if (!checkHash())
		return 0;

	return 1;
}

int XMLQtInterface::isValid() const
{
	QXmlSchema schema;
	if(_schemaName.length() > 0)
		schema.load( QUrl::fromLocalFile((QString::fromStdString(_schemaName))) );

	if ( schema.isValid() )
	{
		QXmlSchemaValidator validator( schema );
		if ( validator.validate( _fileData ) )
			return 1;
		else
		{
			INFO("XMLQtInterface::isValid(): XML file %s is invalid (in reference to schema %s).",
			     _fileName.toStdString().c_str(), _schemaName.c_str());
		}
	}
	else
	{
		QXmlSchemaValidator validator;
		if ( validator.validate( _fileData ) )
			return 1;
		else
		{
			INFO("XMLQtInterface::isValid(): XML file %s is invalid (in reference to its schema).",
			     _fileName.toStdString().c_str());
		}
	}
	return 0;
}

int XMLQtInterface::insertStyleFileDefinition(const QString &fileName) const
{
	std::string path = fileName.toStdString();
	std::fstream stream(path.c_str());
	std::string line;
	std::string styleDef("\n<?xml-stylesheet type=\"text/xsl\" href=\"OpenGeoSysGLI.xsl\"?>");

	if (!stream.is_open())
	{
		WARN("XMLQtInterface::insertStyleFileDefinition(): Could not open file %s.",
		     path.c_str());
		return 0;
	}

	stream.seekp(43 * sizeof(char),std::ios_base::beg); // go to the correct position in the stream
	stream.write(styleDef.c_str(), 60 * sizeof(char)); // write new line with xml-stylesheet definition
	stream.close();
	return 1;
}

bool XMLQtInterface::checkHash() const
{
	QString md5FileName(_fileName + ".md5");
	QByteArray fileHash = QCryptographicHash::hash(_fileData, QCryptographicHash::Md5);

	QFile file(md5FileName);
	if (file.open(QIODevice::ReadOnly))
	{
		if(file.readAll() == fileHash)
			return true;
		INFO("Hashfile does not match data ... checking file ...");
	}

	if (!this->isValid())
		return false;

	QFile fileMD5(md5FileName);
	if(fileMD5.open(QIODevice::WriteOnly))
	{
		fileMD5.write(fileHash);
		fileMD5.close();
		INFO("File is valid, hashfile written.");
	}
	else
		WARN("File is valid but could not write hashfile!");
	return true;
}

bool XMLQtInterface::isHashGood(const QByteArray &hash) const
{
	QByteArray fileHash = QCryptographicHash::hash(_fileData, QCryptographicHash::Md5);
	if(hash != fileHash)
	{
		INFO("Hashfile does not match data ... checking file ...");
		return false;
	}
	return true;
}

} // end namespace FileIO
