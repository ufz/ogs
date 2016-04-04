/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "XmlPrjInterface.h"

#include <ostream>
#include <vector>

#include <logog/include/logog.hpp>

#include "Applications/ApplicationsLib/ProjectData.h"
#include "GeoLib/GEOObjects.h"

#include "BaseLib/FileTools.h"
#include "BaseLib/FileFinder.h"

#include "FileIO/Writer.h"
#include "FileIO/readMeshFromFile.h"
#include "FileIO/VtkIO/VtuInterface.h"
#include "FileIO/XmlIO/Qt/XmlGmlInterface.h"
#include "FileIO/XmlIO/Qt/XmlStnInterface.h"

#include "MeshLib/Mesh.h"

#include <QFile>
#include <QFileInfo>
#include <QtXml/QDomDocument>

namespace FileIO
{
XmlPrjInterface::XmlPrjInterface(ProjectData& project) :
	XMLInterface(), XMLQtInterface(BaseLib::FileFinder().getPath("OpenGeoSysProject.xsd")),
	_project(project)
{
}

int XmlPrjInterface::readFile(const QString &fileName)
{
	if (XMLQtInterface::readFile(fileName) == 0)
		return 0;

	QFileInfo fi(fileName);
	QString const path = (fi.path().length() > 3) ? QString(fi.path() + "/") : fi.path();

	QDomDocument doc("OGS-PROJECT-DOM");
	doc.setContent(_fileData);
	QDomElement const docElement = doc.documentElement(); //OpenGeoSysProject
	if (docElement.nodeName().compare("OpenGeoSysProject"))
	{
		ERR("Unexpected XML root.");
		return 0;
	}

	QDomNodeList const fileList = docElement.childNodes();

	for (int i = 0; i < fileList.count(); ++i)
	{
		QString const file_node(fileList.at(i).nodeName());
		if (file_node.compare("input") == 0)
		{
			if (int n_errors = readInputFiles(fileList.at(i), path))
				ERR ("Error reading %d input files.", n_errors);
		}
	}
	return 0;
}

int XmlPrjInterface::readInputFiles(QDomNode const& node, QString const& path)
{
	QDomNodeList const fileList = node.childNodes();
	int result (0);

	for (int i = 0; i < fileList.count(); i++)
	{
		QString const file_node(fileList.at(i).nodeName());
		if (file_node.compare("geometry") == 0)
		{
			XmlGmlInterface gml(*(_project.getGEOObjects()));
			const QDomNodeList childList = fileList.at(i).childNodes();
			for(int j = 0; j < childList.count(); j++)
			{
				const QDomNode child_node (childList.at(j));
				if (child_node.nodeName().compare("file") == 0)
					result += gml.readFile(QString(path + child_node.toElement().text()));
			}
		}
		else if (file_node.compare("stations") == 0)
		{
			XmlStnInterface stn(*(_project.getGEOObjects()));
			const QDomNodeList childList = fileList.at(i).childNodes();
			for(int j = 0; j < childList.count(); j++)
				if (childList.at(j).nodeName().compare("file") == 0)
					result += stn.readFile(QString(path + childList.at(j).toElement().text()));
		}
		else if (file_node.compare("mesh") == 0)
		{
			const std::string msh_name = path.toStdString() + fileList.at(i).toElement().text().toStdString();
			MeshLib::Mesh* mesh = FileIO::readMeshFromFile(msh_name);
			if (mesh)
				_project.addMesh(mesh);
			else
				result++;
		}
	}

	return result;
}

int XmlPrjInterface::writeToFile(const std::string& filename)
{
	_filename = filename;
	return FileIO::Writer::writeToFile(filename);
}

bool XmlPrjInterface::write()
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
	QDomElement inputTag = doc.createElement("input");
	root.appendChild(inputTag);

	// MSH
	std::vector<MeshLib::Mesh*> const& mesh_vec = _project.getMeshObjects();
	for (MeshLib::Mesh* mesh : mesh_vec)
	{
		// write mesh file
		std::string const name(path + mesh->getName() + ".vtu");
		FileIO::VtuInterface vtkIO(mesh, 0, false);
		if (vtkIO.writeToFile(name.c_str()))
		{
			// write entry in project file
			QDomElement meshTag = doc.createElement("mesh");
			inputTag.appendChild(meshTag);
			QDomElement fileNameTag = doc.createElement("file");
			meshTag.appendChild(fileNameTag);
			QDomText const fileNameText =
				doc.createTextNode(QString::fromStdString(mesh->getName() + ".vtu"));
			fileNameTag.appendChild(fileNameText);
		}
		else
			ERR("Error writing mesh-file \"%s\".", name.c_str());
	}

	// GML
	std::vector<std::string> geo_names;
	geoObjects->getGeometryNames(geo_names);
	for (std::string name : geo_names)
	{
		// write GLI file
		XmlGmlInterface gml(*geoObjects);
		gml.setNameForExport(name);
		if (gml.writeToFile(std::string(path + name + ".gml")) == 0)
		{
			// write entry in project file
			QDomElement geoTag = doc.createElement("geometry");
			inputTag.appendChild(geoTag);
			QDomElement fileNameTag = doc.createElement("file");
			geoTag.appendChild(fileNameTag);
			QDomText const fileNameText = doc.createTextNode(QString::fromStdString(name + ".gml"));
			fileNameTag.appendChild(fileNameText);
		}
		else
			ERR("Error writing geometry-file \"%s\".", name.c_str());
	}

	// STN
	std::vector<std::string> stn_names;
	geoObjects->getStationVectorNames(stn_names);
	for (std::string name : stn_names)
	{
		// write STN file
		XmlStnInterface stn(*geoObjects);
		stn.setNameForExport(name);

		if (stn.writeToFile(path + name + ".stn"))
		{
			// write entry in project file
			QDomElement geoTag = doc.createElement("stations");
			inputTag.appendChild(geoTag);
			QDomElement fileNameTag = doc.createElement("file");
			geoTag.appendChild(fileNameTag);
			QDomText const fileNameText = doc.createTextNode(QString::fromStdString(name + ".stn"));
			fileNameTag.appendChild(fileNameText);
		}
		else
			ERR("XmlPrjInterface::writeFile(): Error writing stations-file \"%s\".", name.c_str());
	}

	std::string xml = doc.toString().toStdString();
	_out << xml;
	return true;
}
}
