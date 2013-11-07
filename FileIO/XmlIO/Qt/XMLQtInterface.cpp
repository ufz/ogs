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
#include <QCryptographicHash>

// ThirdParty/logog
#include "logog/include/logog.hpp"

namespace FileIO
{

XMLQtInterface::XMLQtInterface(const std::string &schemaFile) :
		_schemaName(schemaFile)
{}

int XMLQtInterface::isValid(const QString &fileName) const
{
	QXmlSchema schema;
	if(_schemaName.length() > 0)
		schema.load( QUrl::fromLocalFile((QString::fromStdString(_schemaName))) );

	if ( schema.isValid() )
	{
		QXmlSchemaValidator validator( schema );
		if ( validator.validate( QUrl::fromLocalFile((fileName))) )
			return 1;
		else
		{
			INFO("XMLQtInterface::isValid(): XML file %s is invalid (in reference to schema %s).",
			     fileName.data(), _schemaName.c_str());
			return 0;
		}
	}
	else
	{
		QXmlSchemaValidator validator;
		if ( validator.validate( QUrl::fromLocalFile((fileName))) )
			return 1;
		else
		{
			INFO("XMLQtInterface::isValid(): XML file %s is invalid (in reference to its schema).",
			     fileName.data(), _schemaName.c_str());
			return 0;
		}
	}
}

void XMLQtInterface::setSchema(const std::string &schemaName)
{
	_schemaName = schemaName;
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

bool XMLQtInterface::checkHash(const QString &fileName) const
{
	QFileInfo fi(fileName);
	QString md5FileName(fileName + ".md5");

	std::ifstream md5( md5FileName.toStdString().c_str() );
	if (md5.is_open())
	{
		char* md5HashStr = new char[16];
		md5.read(md5HashStr, 16);
		QByteArray md5Hash(md5HashStr, 16);
		delete[] md5HashStr;
		if (isHashGood(fileName, md5Hash))
			return true;
	}

	if (!this->isValid(fileName))
		return false;

	INFO("File is valid, writing hashfile.");
	QByteArray hash = calcHash(fileName);
	std::ofstream out( md5FileName.toStdString().c_str(), std::ios::out );
	out.write(hash.data(), 16);
	out.close();
	return true;
}

bool XMLQtInterface::isHashGood(const QString &fileName, const QByteArray &hash) const
{
	int hashLength = hash.length();
	QByteArray fileHash = calcHash(fileName);
	if (fileHash.length() != hashLength)
		return false;
	for (int i = 0; i < hashLength; i++)
		if (fileHash[i] != hash[i])
		{
			INFO("Hashfile does not match data ... checking file ...");
			return false;
		}
	return true;
}

QByteArray XMLQtInterface::calcHash(const QString &fileName) const
{
	std::ifstream is(fileName.toStdString().c_str(), std::ios::binary );
	is.seekg (0, std::ios::end);
	int length = is.tellg();
	is.seekg (0, std::ios::beg);
	char* buffer = new char [length];
	is.read (buffer,length);
	is.close();

	QByteArray hash = QCryptographicHash::hash(buffer, QCryptographicHash::Md5);
	delete [] buffer;
	return hash;
}

} // end namespace FileIO
