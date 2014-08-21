/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Implementation of the XmlProjectInterface class.
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

#include "XmlProjectInterface.h"

#include "XmlGmlInterface.h"
#include "XmlStnInterface.h"
#include "XmlCndInterface.h"
#include "XmlIO/Boost/BoostVtuInterface.h"

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
XmlProjectInterface::XmlProjectInterface(ProjectData& project) :
	XMLInterface(), XMLQtInterface(BaseLib::FileFinder().getPath("OpenGeoSysProject.xsd")),
	_project(project)
{
}

int XmlProjectInterface::readFile(const QString &fileName)
{
	if(XMLQtInterface::readFile(fileName) == 0)
		return 0;

	QFileInfo const fi(fileName);
	QString const path = (fi.path().length() > 3) ? QString(fi.path() + "/") : fi.path();

	QDomDocument doc("OGS-PROJECT-DOM");
	doc.setContent(_fileData);
	QDomElement const docElement = doc.documentElement(); //OpenGeoSysProject
	if (docElement.nodeName().compare("OpenGeoSysProject"))
	{
		ERR("XmlProjectInterface::readFile(): Unexpected XML root.");
		return 0;
	}

	QDomElement file_node = docElement.firstChildElement();
	while (!file_node.isNull())
	{
		if (file_node.nodeName().compare("geofile") == 0)
		{
			XmlGmlInterface gml(*(_project.getGEOObjects()));
			gml.readFile(QString(path + file_node.toElement().text()));
		}
		else if (file_node.nodeName().compare("stnfile") == 0)
		{
			XmlStnInterface stn(*(_project.getGEOObjects()));
			stn.readFile(QString(path + file_node.toElement().text()));
		}
		else if (file_node.nodeName().compare("mshfile") == 0)
		{
			std::string const msh_name = 
				path.toStdString() + file_node.toElement().text().toStdString();
			MeshLib::Mesh* mesh = FileIO::readMeshFromFile(msh_name);
			if (mesh)
				_project.addMesh(mesh);
		}
		else if (file_node.nodeName().compare("process") == 0)
		{
			QDomElement process_node = file_node.firstChildElement();
			while (!process_node.isNull())
			{
				if (process_node.nodeName().compare("cndfile") == 0)
				{
					const std::string cnd_name = 
						path.toStdString() + process_node.toElement().text().toStdString();
					XmlCndInterface cnd(_project);
					cnd.readFile(cnd_name);
				}
				if (file_node.nodeName().compare("matfile") == 0)
				{
					INFO ("Material file reader not yet implemented.");
				}
				if (file_node.nodeName().compare("numfile") == 0)
				{
					INFO ("Numeric file reader not yet implemented.");
				}
				process_node = process_node.nextSiblingElement();
			}
		}
		file_node = file_node.nextSiblingElement();
	}

	return 1;
}

int XmlProjectInterface::writeToFile(std::string filename)
{
	_filename = filename;
	return FileIO::Writer::writeToFile(filename);
}

bool XmlProjectInterface::write()
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

	
	{ // GML
	std::vector<std::string> geoNames;
	geoObjects->getGeometryNames(geoNames);
	XmlGmlInterface gml(*geoObjects);
	for (std::vector<std::string>::const_iterator it(geoNames.begin()); it != geoNames.end(); ++it)
	{
		gml.setNameForExport(*it);
		this->writeFile(doc, root, "geofile", gml, path, *it, "gml");
	}
	} // end GML-scope

	{ // MSH
	std::vector<MeshLib::Mesh*> const msh_vec = _project.getMeshObjects();
	BoostVtuInterface msh;
	for (std::vector<MeshLib::Mesh*>::const_iterator it(msh_vec.begin()); it != msh_vec.end(); ++it)
	{
		std::string const name((*it)->getName());
		msh.setMesh(*it);
		if (msh.writeToFile(path + name + ".vtu"))
		{
			// write entry in project file
			QDomElement fileTag = doc.createElement("mshfile");
			root.appendChild(fileTag);
			QDomText const fileNameText =
				doc.createTextNode(QString::fromStdString(name + ".vtu"));
			fileTag.appendChild(fileNameText);
		}
		else
			ERR("XmlProjectInterface::writeFile(): Error writing vtu-file \"%s\".", name.c_str());
	}
	} // end MSH-scope

	{ // STN
	std::vector<std::string> stnNames;
	geoObjects->getStationVectorNames(stnNames);
	XmlStnInterface stn(*geoObjects);
	for (std::vector<std::string>::const_iterator it(stnNames.begin()); it != stnNames.end(); ++it)
	{
		stn.setNameForExport(*it);
		this->writeFile(doc, root, "stnfile", stn, path, *it, "stn");
	}
	} // end STN-scope

	// process-specific data
	QDomElement processTag = doc.createElement("process");
	std::string const process_name ("Undefined"); // TODO: extract process name from project data
	processTag.setAttribute( "type", QString::fromStdString(process_name) );
	root.appendChild(processTag);

	{ // CND
	const std::vector<FEMCondition*> &cnd_vec (_project.getConditions());
	XmlCndInterface cnd(_project);
	if (!cnd_vec.empty())
	{		
		std::string const cnd_name (process_name + ".cnd");
		this->writeFile(doc, processTag, "cndfile", cnd, path, process_name, "cnd");
	}
	} // end CND-scope

	// TODO: writer for MAT-files

	// TODO: writer for NUM-files

	// if no process-dependent information exists, the process tag can be removed again
	if (processTag.toElement().childNodes().isEmpty())
		root.removeAttribute("process");

	std::string xml = doc.toString().toStdString();
	_out << xml;
	return true;
}


void XmlProjectInterface::writeFile(QDomDocument &doc,
                                    QDomElement &parent,
                                    QString const& tag_name,
                                    XMLInterface &xml_interface,
                                    std::string const& path,
                                    std::string const& file_name,
                                    std::string const& extension)
{
	if (xml_interface.writeToFile(path + file_name + "." + extension))
	{
		// write entry in project file
		QDomElement fileTag = doc.createElement(tag_name);
		parent.appendChild(fileTag);
		QDomText const fileNameText =
			doc.createTextNode(QString::fromStdString(file_name + "." + extension));
		fileTag.appendChild(fileNameText);
	}
	else
		ERR("XmlProjectInterface::writeFile(): Error writing %s-file \"%s\".", extension.c_str(), file_name.c_str());
}


}
