/**
 * \file XMLInterface.cpp
 * 18/02/2010 KR Initial implementation
 */

#include "ProjectData.h"
#include "XMLInterface.h"

//#include "Configure.h"

#include <QCryptographicHash>
#include <QFileInfo>
#ifdef QT_USE_QTXMLPATTERNS
#include <QtXmlPatterns/QXmlSchema>
#include <QtXmlPatterns/QXmlSchemaValidator>
#endif // QT_USE_QTXMLPATTERNS

namespace FileIO
{

XMLInterface::XMLInterface(ProjectData* project, const std::string &schemaFile)
: _project(project), _exportName(""), _schemaName(schemaFile)
{
}

int XMLInterface::isValid(const QString &fileName) const
{
#ifdef QT_USE_QTXMLPATTERNS
	QXmlSchema schema;
	schema.load( QUrl::fromLocalFile((QString::fromStdString(_schemaName))) );

	if ( schema.isValid() )
	{
		QXmlSchemaValidator validator( schema );
		if ( validator.validate( QUrl::fromLocalFile((fileName))) )
			return 1;
		else
		{
			std::cout <<
			"XMLInterface::isValid() - XML File is invalid (in reference to schema " <<
			_schemaName << ")." << std::endl;
			return 0;
		}
	}
	else
	{
		std::cout << "XMLInterface::isValid() - Schema " << _schemaName <<
		" is invalid." << std::endl;
		return 0;
	}
#else // ifdef QT_USE_QTXMLPATTERNS
	Q_UNUSED (fileName);
	std::cout <<
	"XMLInterface: XML schema validation skipped. Qt 4.6 is required for validation." <<
	std::endl;
	return 1;
#endif // QT_USE_QTXMLPATTERNS
}

void XMLInterface::setSchema(const std::string &schemaName)
{
	_schemaName = schemaName;
}


int XMLInterface::insertStyleFileDefinition(const QString &fileName) const
{
	std::string path = fileName.toStdString();
	std::fstream stream(path.c_str());
	std::string line;
	std::string styleDef("\n<?xml-stylesheet type=\"text/xsl\" href=\"OpenGeoSysGLI.xsl\"?>");

	if (!stream.is_open())
	{
		std::cout << "XMLInterface::insertStyleFileDefinition() - Could not open file...\n";
		return 0;
	}

	stream.seekp(43 * sizeof(char),std::ios_base::beg); // go to the correct position in the stream
	stream.write(styleDef.c_str(), 60 * sizeof(char)); // write new line with xml-stylesheet definition
	stream.close();
	return 1;
}

bool XMLInterface::checkHash(const QString &fileName) const
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
		if (hashIsGood(fileName, md5Hash))
			return true;
	}

	if (!this->isValid(fileName))
		return false;

	std::cout << "File is valid, writing hashfile..." << std::endl;
	QByteArray hash = calcHash(fileName);
	std::ofstream out( md5FileName.toStdString().c_str(), std::ios::out );
	out.write(hash.data(), 16);
	out.close();
	return true;
}

bool XMLInterface::hashIsGood(const QString &fileName, const QByteArray &hash) const
{
	int hashLength = hash.length();
	QByteArray fileHash = calcHash(fileName);
	if (fileHash.length() != hashLength)
		return false;
	for (int i = 0; i < hashLength; i++)
		if (fileHash[i] != hash[i])
		{
			std::cout << "Hashfile does not match data ... checking file ..." <<
			std::endl;
			return false;
		}
	return true;
}

QByteArray XMLInterface::calcHash(const QString &fileName) const
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

}
