/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Implementation of the XmlGspInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "XmlGspInterface.h"

#include "XmlCndInterface.h"
#include "XmlGmlInterface.h"
#include "XmlStnInterface.h"

#include "FileTools.h"
#include "FileFinder.h"
#include "FileIO/Legacy/MeshIO.h"
#include "FileIO/readMeshFromFile.h"
#include "Mesh.h"

#include <QFile>
#include <QFileInfo>
#include <QtXml/QDomDocument>

namespace FileIO
{
XmlGspInterface::XmlGspInterface(ProjectData& project) :
	XMLInterface(), XMLQtInterface(BaseLib::FileFinder().getPath("OpenGeoSysProject.xsd")),
	_project(project)
{
}

int XmlGspInterface::readFile(const QString &fileName)
{
	if(XMLQtInterface::readFile(fileName) == 0)
		return 0;

	QFileInfo fi(fileName);
	QString path = (fi.path().length() > 3) ? QString(fi.path() + "/") : fi.path();

	QDomDocument doc("OGS-PROJECT-DOM");
	doc.setContent(_fileData);
	QDomElement docElement = doc.documentElement(); //OpenGeoSysProject
	if (docElement.nodeName().compare("OpenGeoSysProject"))
	{
		ERR("XmlGspInterface::readFile(): Unexpected XML root.");
		return 0;
	}

	QDomNodeList fileList = docElement.childNodes();

	for(int i = 0; i < fileList.count(); i++)
	{
		const QString file_node(fileList.at(i).nodeName());
		if (file_node.compare("geo") == 0)
		{
			XmlGmlInterface gml(*(_project.getGEOObjects()));
			const QDomNodeList childList = fileList.at(i).childNodes();
			for(int j = 0; j < childList.count(); j++)
			{
				const QDomNode child_node (childList.at(j));
				if (child_node.nodeName().compare("file") == 0)
				{
					DBUG("XmlGspInterface::readFile(): path: \"%s\".",
					     path.data());
					DBUG("XmlGspInterface::readFile(): file name: \"%s\".",
					     (child_node.toElement().text()).data());
					gml.readFile(QString(path + child_node.toElement().text()));
				}
			}
		}
		else if (file_node.compare("stn") == 0)
		{
			XmlStnInterface stn(*(_project.getGEOObjects()));
			const QDomNodeList childList = fileList.at(i).childNodes();
			for(int j = 0; j < childList.count(); j++)
				if (childList.at(j).nodeName().compare("file") == 0)
					stn.readFile(QString(path +
					                     childList.at(j).toElement().text()));
		}
		else if (file_node.compare("msh") == 0)
		{
			const std::string msh_name = path.toStdString() +
			                             fileList.at(i).toElement().text().toStdString();
			MeshLib::Mesh* msh = FileIO::readMeshFromFile(msh_name);
			if (msh)
				_project.addMesh(msh);
		}
		else if (file_node.compare("cnd") == 0)
		{
			const std::string cnd_name = path.toStdString() +
			                             fileList.at(i).toElement().text().toStdString();
			XmlCndInterface cnd(_project);
			cnd.readFile(cnd_name);
		}
	}

	return 1;
}

int XmlGspInterface::writeToFile(std::string filename)
{
	_filename = filename;
	return FileIO::Writer::writeToFile(filename);
}

bool XmlGspInterface::write()
{
	GeoLib::GEOObjects* geoObjects = _project.getGEOObjects();
	QFileInfo fi(QString::fromStdString(_filename));
	std::string path((fi.absolutePath()).toStdString() + "/");

	_out << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"; // xml definition
	_out << "<?xml-stylesheet type=\"text/xsl\" href=\"OpenGeoSysProject.xsl\"?>\n\n"; // stylefile definition

	QDomDocument doc("OGS-PROJECT-DOM");
	QDomElement root = doc.createElement("OpenGeoSysProject");
	root.setAttribute( "xmlns:ogs", "http://www.opengeosys.org" );
	root.setAttribute( "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance" );
	root.setAttribute( "xsi:noNamespaceSchemaLocation",
	                   "http://www.opengeosys.org/images/xsd/OpenGeoSysProject.xsd" );

	doc.appendChild(root);

	// GML
	std::vector<std::string> geoNames;
	geoObjects->getGeometryNames(geoNames);
	for (std::vector<std::string>::const_iterator it(geoNames.begin()); it != geoNames.end(); ++it)
	{
		// write GLI file
		XmlGmlInterface gml(*geoObjects);
		std::string name(*it);
		gml.setNameForExport(name);
		if (gml.writeToFile(std::string(path + name + ".gml")))
		{
			// write entry in project file
			QDomElement geoTag = doc.createElement("geo");
			root.appendChild(geoTag);
			QDomElement fileNameTag = doc.createElement("file");
			geoTag.appendChild(fileNameTag);
			QDomText fileNameText =
			        doc.createTextNode(QString::fromStdString(name + ".gml"));
			fileNameTag.appendChild(fileNameText);
		}
	}

	// MSH
	const std::vector<MeshLib::Mesh*> msh_vec = _project.getMeshObjects();
	for (std::vector<MeshLib::Mesh*>::const_iterator it(msh_vec.begin()); it != msh_vec.end(); ++it)
	{
		// write mesh file
		Legacy::MeshIO meshIO;
		meshIO.setMesh(*it);
		std::string fileName(path + (*it)->getName());
		meshIO.writeToFile(fileName);

		// write entry in project file
		QDomElement mshTag = doc.createElement("msh");
		root.appendChild(mshTag);
		QDomElement fileNameTag = doc.createElement("file");
		mshTag.appendChild(fileNameTag);
		QDomText fileNameText = doc.createTextNode(QString::fromStdString((*it)->getName()));
		fileNameTag.appendChild(fileNameText);
	}

	// STN
	std::vector<std::string> stnNames;
	geoObjects->getStationVectorNames(stnNames);
	for (std::vector<std::string>::const_iterator it(stnNames.begin()); it != stnNames.end(); ++it)
	{
		// write STN file
		XmlStnInterface stn(*geoObjects);
		std::string name(*it);
		stn.setNameForExport(name);

		if (stn.writeToFile(path + name + ".stn"))
		{
			// write entry in project file
			QDomElement geoTag = doc.createElement("stn");
			root.appendChild(geoTag);
			QDomElement fileNameTag = doc.createElement("file");
			geoTag.appendChild(fileNameTag);
			QDomText fileNameText =
			        doc.createTextNode(QString::fromStdString(name + ".stn"));
			fileNameTag.appendChild(fileNameText);
		}
		else
			ERR("XmlGspInterface::writeFile(): Error writing stn-file \"%s\".", name.c_str());
	}

	// CND
	const std::vector<FEMCondition*> &cnd_vec (_project.getConditions());
	if (!cnd_vec.empty())
	{
		XmlCndInterface cnd(_project);
		const std::string cnd_name (BaseLib::extractBaseNameWithoutExtension(_filename) + ".cnd");
		if (cnd.writeToFile(path + cnd_name))
		{
			// write entry in project file
			QDomElement cndTag = doc.createElement("cnd");
			root.appendChild(cndTag);
			QDomElement fileNameTag = doc.createElement("file");
			cndTag.appendChild(fileNameTag);
			QDomText fileNameText = doc.createTextNode(QString::fromStdString(cnd_name));
			fileNameTag.appendChild(fileNameText);
		}
		else
			ERR("XmlGspInterface::writeFile(): Error writing cnd-file \"%s\".", cnd_name.c_str());
	}

	std::string xml = doc.toString().toStdString();
	_out << xml;
	return true;
}
}
