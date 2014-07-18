/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-23
 * \brief  Implementation of the XmlCndInterface class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "XmlCndInterface.h"

#include <QFile>
#include <QTextCodec>
#include <QtXml/QDomDocument>
#include <QStringList>

#include "FEMCondition.h"
#include "ProjectData.h"
#include "FileFinder.h"


namespace FileIO
{
XmlCndInterface::XmlCndInterface(ProjectData &project)
	: XMLInterface(), XMLQtInterface(BaseLib::FileFinder().getPath("OpenGeoSysCND.xsd")),
	  _type(FEMCondition::UNSPECIFIED), _project(project)
{
}

int XmlCndInterface::readFile(const QString &fileName)
{
	if(XMLQtInterface::readFile(fileName) == 0)
		return 0;

	QDomDocument doc("OGS-Cond-DOM");
	doc.setContent(_fileData);
	QDomElement docElement = doc.documentElement(); //root element, used for identifying file-type
	if (docElement.nodeName().compare("OpenGeoSysCond"))
	{
		ERR("XMLInterface::readFEMCondFile(): Unexpected XML root.");
		return 0;
	}

	std::size_t const n_cond_before(this->_project.getConditions().size());
	QDomNodeList lists = docElement.childNodes();
	for (int i = 0; i < lists.count(); i++)
	{
		const QDomNode list_node (lists.at(i));
		if (list_node.nodeName().compare("BoundaryConditions") == 0)
			readConditions(list_node, FEMCondition::BOUNDARY_CONDITION);
		else if (list_node.nodeName().compare("InitialConditions") == 0)
			readConditions(list_node, FEMCondition::INITIAL_CONDITION);
		else if (list_node.nodeName().compare("SourceTerms") == 0)
			readConditions(list_node, FEMCondition::SOURCE_TERM);
	}
	std::size_t const n_cond_after(this->_project.getConditions().size());
	if (n_cond_after-n_cond_before > 0)
		return 1;     //do something like _geoObjects->addStationVec(stations, stnName, color);
	else
	{
		WARN("XMLInterface::readFEMCondFile(): No FEM Conditions found.");
		return 0;
	}

	return 1;
}

void XmlCndInterface::readConditions(const QDomNode &listRoot,
                                     FEMCondition::CondType type)
{
	QDomElement cond = listRoot.firstChildElement();
	while (!cond.isNull())
	{
		std::string geometry_name ( cond.attribute("geometry").toStdString() );
		if (this->_project.getGEOObjects()->exists(geometry_name) >= 0 ||
		    this->_project.meshExists(geometry_name))
		{
			FEMCondition* c ( new FEMCondition(geometry_name, type) );

			QDomNodeList condProperties = cond.childNodes();
			for (int i = 0; i < condProperties.count(); i++)
			{
				const QDomNode prop_node (condProperties.at(i));
				if (condProperties.at(i).nodeName().compare("Process") == 0)
				{
					QDomNodeList processProps = prop_node.childNodes();
					for (int j = 0; j < processProps.count(); j++)
					{
						const QString prop_name(processProps.at(j).nodeName());
						if (prop_name.compare("Type") == 0)
							c->setProcessType(FiniteElement::convertProcessType(processProps.at(j).toElement().text().toStdString()));
						else if (prop_name.compare("Variable") == 0)
							c->setProcessPrimaryVariable(FiniteElement::convertPrimaryVariable(processProps.at(j).toElement().text().toStdString()));
					}
				}
				else if (prop_node.nodeName().compare("Geometry") == 0)
				{
					QDomNodeList geoProps = prop_node.childNodes();
					GeoLib::GEOTYPE geo_type;
					std::string geo_obj_name;
					for (int j = 0; j < geoProps.count(); j++)
					{
						const QString prop_name(geoProps.at(j).nodeName());
						if (prop_name.compare("Type") == 0)
							geo_type = GeoLib::convertGeoType(geoProps.at(j).toElement().text().toStdString());
						else if (prop_name.compare("Name") == 0)
							geo_obj_name = geoProps.at(j).toElement().text().toStdString();
					}
					c->initGeometricAttributes(geometry_name, geo_type, geo_obj_name, *(_project.getGEOObjects()));
				}
				else if (prop_node.nodeName().compare("Distribution") == 0)
				{
					QDomNodeList distProps = prop_node.childNodes();
					for (int j = 0; j < distProps.count(); j++)
					{
						const QString prop_name(distProps.at(j).nodeName());
						if (prop_name.compare("Type") == 0)
							c->setProcessDistributionType(FiniteElement::convertDisType(distProps.at(j).toElement().text().toStdString()));
						else if (prop_name.compare("Value") == 0)
						{
							std::vector<size_t> disNodes;
							std::vector<double> disValues;
							if (c->getProcessDistributionType()==FiniteElement::CONSTANT ||
								c->getProcessDistributionType()==FiniteElement::CONSTANT_NEUMANN ||
								c->getProcessDistributionType()==FiniteElement::NODESCONSTANT)
								disValues.push_back( distProps.at(j).toElement().text().toDouble() );
							else if (c->getProcessDistributionType()==FiniteElement::LINEAR ||
								     c->getProcessDistributionType()==FiniteElement::LINEAR_NEUMANN ||
									 c->getProcessDistributionType()==FiniteElement::DIRECT)
							{
								QString text = distProps.at(j).toElement().text();
								QStringList list = text.split(QRegExp("\\t"));
								size_t count(0);
								for (QStringList::iterator it=list.begin(); it!=list.end(); ++it)
								{
									QString val (it->trimmed());
									if (!val.isEmpty())
									{
										if (count%2==0)
											disNodes.push_back(val.toInt());
										else
											disValues.push_back(val.toDouble());
										count++;
									}
								}
							}
							else
								ERR("XmlCndInterface::readConditions(): Distribution type not supported.");
							c->setDisValues(disNodes, disValues);
						}
					}
				}
			}
			this->_project.addCondition(c);
		}
		else
		{
			ERR("XmlCndInterface::readConditions(): No geometry \"%s\" found.", geometry_name.c_str());
		}
		cond = cond.nextSiblingElement();
	}
}

bool XmlCndInterface::write()
{
	_out << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"; // xml definition
	_out << "<?xml-stylesheet type=\"text/xsl\" href=\"OpenGeoSysCND.xsl\"?>\n\n"; // stylefile definition

	QDomDocument doc("OGS-CND-DOM");
	QDomElement root = doc.createElement("OpenGeoSysCond");
	root.setAttribute( "xmlns:ogs", "http://www.opengeosys.org" );
	root.setAttribute( "xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance" );
	root.setAttribute( "xsi:noNamespaceSchemaLocation", "http://www.opengeosys.org/images/xsd/OpenGeoSysCND.xsd" );

	std::vector<FEMCondition*> const& conditions (_project.getConditions(
	                                                     FiniteElement::INVALID_PROCESS,
	                                                     _exportName,
	                                                     _type) );

	if (conditions.empty())
		return 1;

	const QString export_name (QString::fromStdString(_exportName));
	doc.appendChild(root);

	size_t nConditions (conditions.size());

	// create root nodes for various conditions types if there is at least one condition of that type
	QDomElement ic_root, bc_root, st_root;
	if (std::count_if(conditions.begin(), conditions.end(), 
		[](FEMCondition const*const cond){return cond->getCondType() == FEMCondition::INITIAL_CONDITION;}))
		ic_root = this->getCondListElement(doc, root, "InitialConditions");
	if (std::count_if(conditions.begin(), conditions.end(), 
		[](FEMCondition const*const cond){return cond->getCondType() == FEMCondition::BOUNDARY_CONDITION;}))
		bc_root = this->getCondListElement(doc, root, "BoundaryConditions");
	if (std::count_if(conditions.begin(), conditions.end(), 
		[](FEMCondition const*const cond){return cond->getCondType() == FEMCondition::SOURCE_TERM;}))
		st_root = this->getCondListElement(doc, root, "SourceTerms");

	for (size_t i = 0; i < nConditions; i++)
	{
		FEMCondition::CondType current_type = conditions[i]->getCondType();
		if (current_type == _type || _type == FEMCondition::UNSPECIFIED)
		{
			QDomElement listTag;
			QString condText;

			if (current_type == FEMCondition::BOUNDARY_CONDITION)
			{
				listTag = bc_root;
				condText = "BC";
			}
			else if (current_type == FEMCondition::INITIAL_CONDITION)
			{
				listTag = ic_root;
				condText = "IC";
			}
			else if (current_type == FEMCondition::SOURCE_TERM)
			{
				listTag = st_root;
				condText = "ST";
			}
			else
			{
				ERR("XmlCndInterface::writeFile(): Unspecified FEMConditions found ... Abort writing.");
				return 0;
			}
			this->writeCondition(doc, listTag, conditions[i], condText, export_name);
		}
	}
	std::string xml = doc.toString().toStdString();
	_out << xml;

	return true;
}

void XmlCndInterface::writeCondition(QDomDocument doc, QDomElement &listTag,
                                     const FEMCondition* cond, const QString &condText,
                                     const QString &geometryName) const
{
	QString geoName (QString::fromStdString(cond->getAssociatedGeometryName()));

	if ((geometryName.length() > 0) && (geoName.compare(geometryName) != 0))
	{
		WARN("XmlCndInterface::writeCondition(): Geometry name not matching, skipping condition \"%s\".",
		     cond->getGeoName().c_str());
		return;
	}

	QDomElement condTag ( doc.createElement(condText) );
	condTag.setAttribute("geometry", geoName);
	listTag.appendChild(condTag);

	QDomElement processTag ( doc.createElement("Process") );
	condTag.appendChild(processTag);
	QDomElement processTypeTag ( doc.createElement("Type") );
	processTag.appendChild(processTypeTag);
	QDomText processTypeText ( doc.createTextNode(
		QString::fromStdString(FiniteElement::convertProcessTypeToString(cond->getProcessType()))) );
	processTypeTag.appendChild(processTypeText);
	QDomElement processPVTag ( doc.createElement("Variable") );
	processTag.appendChild(processPVTag);
	QDomText processPVText ( doc.createTextNode(
		QString::fromStdString(FiniteElement::convertPrimaryVariableToString(cond->getProcessPrimaryVariable()))) );
	processPVTag.appendChild(processPVText);

	QDomElement geoTag ( doc.createElement("Geometry") );
	condTag.appendChild(geoTag);
	QDomElement geoTypeTag ( doc.createElement("Type") );
	geoTag.appendChild(geoTypeTag);
	QDomText geoTypeText ( doc.createTextNode(
		QString::fromStdString(GeoLib::convertGeoTypeToString(cond->getGeomType()))) );
	geoTypeTag.appendChild(geoTypeText);
	QDomElement geoNameTag ( doc.createElement("Name") );
	geoTag.appendChild(geoNameTag);
	QString geo_obj_name ( QString::fromStdString(cond->getGeoName()) );
	QDomText geoNameText ( doc.createTextNode(geo_obj_name) );
	geoNameTag.appendChild(geoNameText);

	QDomElement disTag ( doc.createElement("Distribution") );
	condTag.appendChild(disTag);
	QDomElement disTypeTag ( doc.createElement("Type") );
	disTag.appendChild(disTypeTag);
	QDomText disTypeText ( doc.createTextNode(
		QString::fromStdString(FiniteElement::convertDisTypeToString(cond->getProcessDistributionType()))) );
	disTypeTag.appendChild(disTypeText);
	QDomElement disValueTag ( doc.createElement("Value") );
	disTag.appendChild(disValueTag);
	/*
	   if (cond->getProcessDistributionType() != FiniteElement::DIRECT)
	   {
	    double dis_value (cond->getDisValue()[0]); //TODO: do this correctly!
	    disValueText = doc.createTextNode(QString::number(dis_value));
	   }
	   else
	    disValueText = doc.createTextNode(QString::fromStdString(cond->getDirectFileName()));
	 */
	const std::vector<size_t> dis_nodes = cond->getDisNodes();
	const std::vector<double> dis_values = cond->getDisValues();
	const size_t nNodes = dis_nodes.size();
	const size_t nValues = dis_values.size();
	std::stringstream ss;
	if (nNodes == 0 && nValues == 1)        // CONSTANT
		ss << dis_values[0];
	else if ((nValues > 0) && (nValues == nNodes)) // LINEAR && DIRECT
	{
		ss << "\n\t";
		for (size_t i = 0; i < nValues; i++)
			ss << dis_nodes[i] << "\t" << dis_values[i] << "\n\t";
	}
	else
	{
		ERR("XmlCndInterface::writeCondition(): Inconsistent length of distribution value array.");
		ss << "-9999";
	}
	std::string dv  = ss.str();
	QDomText disValueText = doc.createTextNode(QString::fromStdString(ss.str()));
	disValueTag.appendChild(disValueText);
}

QDomElement XmlCndInterface::getCondListElement(QDomDocument doc, QDomElement &root,
                                                const QString &text) const
{
	const QDomNodeList list = root.elementsByTagName(text);
	if (list.isEmpty()) {
		QDomElement newListTag(doc.createElement(text));
		root.appendChild(newListTag);
		return newListTag;
	}
	return list.at(0).toElement();
}
}
