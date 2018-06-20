/**
 * \file
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "XmlPrjInterface.h"

#include <iostream>
#include <vector>

#include <QFile>
#include <QFileInfo>
#include <QtXml/QDomDocument>
#include <logog/include/logog.hpp>

#include "BaseLib/BuildInfo.h"
#include "BaseLib/FileFinder.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/IO/Writer.h"

#include "Applications/DataHolderLib/FemCondition.h"

#include "GeoLib/GEOObjects.h"

#include "GeoLib/IO/XmlIO/Qt/XmlGmlInterface.h"
#include "GeoLib/IO/XmlIO/Qt/XmlStnInterface.h"

#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"

namespace FileIO
{
XmlPrjInterface::XmlPrjInterface(DataHolderLib::Project& project)
    : XMLInterface(),
      XMLQtInterface("OpenGeoSysProject.xsd"), _filename(""), _project(project)
{
}

int XmlPrjInterface::readFile(const QString& fileName)
{
    if (XMLQtInterface::readFile(fileName) == 0)
        return 0;

    QFileInfo fi(fileName);
    QString path =
        (fi.path().length() > 3) ? QString(fi.path() + "/") : fi.path();

    QDomNode param_root = QDomNode();
    QDomNode pvar_root = QDomNode();
    QDomDocument doc("OGS-PROJECT-DOM");
    doc.setContent(_fileData);
    QDomElement docElement = doc.documentElement();  // OpenGeoSysProject
    if (docElement.nodeName().compare("OpenGeoSysProject"))
    {
        ERR("XmlPrjInterface::readFile(): Unexpected XML root.");
        return 0;
    }

    QDomNodeList fileList = docElement.childNodes();

    for (int i = 0; i < fileList.count(); i++)
    {
        QDomNode const node(fileList.at(i));
        QString const node_name(node.nodeName());
        QString const file_name(node.toElement().text().trimmed());
        if (file_name.isEmpty())
            continue;

        if (node_name == "geometry")
        {
            GeoLib::IO::XmlGmlInterface gml(_project.getGEOObjects());
            gml.readFile(QString(path + file_name));
        }
        else if (node_name == "stations")
        {
            GeoLib::IO::XmlStnInterface stn(_project.getGEOObjects());
            stn.readFile(QString(path + file_name));
        }
        else if (node_name == "mesh")
        {
            QString const mesh_str(path + file_name);
            std::unique_ptr<MeshLib::Mesh> mesh(
                MeshLib::IO::readMeshFromFile(mesh_str.toStdString()));
            if (mesh)
                _project.addMesh(std::move(mesh));
        }

        else if (node_name == "parameters")
            param_root = node;
        else if (node_name == "process_variables")
            pvar_root = node;
    }

    if (param_root != QDomNode() && pvar_root != QDomNode())
        readConditions(pvar_root, param_root);
    else
        INFO("Skipping process variables");

    return 1;
}

QDomNode XmlPrjInterface::findParam(QDomNode const& param_root,
                                    QString const& param_name) const
{
    QDomNode param = param_root.firstChild();

    while (param != QDomNode())
    {
        QDomNodeList nodeList = param.childNodes();
        for (int i = 0; i < nodeList.count(); i++)
        {
            QDomNode node = nodeList.at(i);
            // std::cout << node.nodeName().toStdString() <<
            // node.toElement().text().toStdString() << std::endl;
            if (node.nodeName() == "name" &&
                node.toElement().text() == param_name)
                return node;
        }
        param = param.nextSibling();
    }
    return QDomNode();
}

void XmlPrjInterface::readConditions(QDomNode const& pvar_root,
                                     QDomNode const& param_root)
{
    QDomNode pvar = pvar_root.firstChild();
    while (pvar != QDomNode())
    {
        DataHolderLib::ProcessVariable process_var;
        QDomNodeList nodeList = pvar.childNodes();
        for (int i = 0; i < nodeList.count(); i++)
        {
            QDomNode const node = nodeList.at(i);
            QString const node_name = node.nodeName();
            if (node_name == "name")
                process_var.name = node.toElement().text().toStdString();
            else if (node_name == "components")
                process_var.components = node.toElement().text().toInt();
            else if (node_name == "order")
                process_var.order = node.toElement().text().toInt();
            else if (node_name == "boundary_conditions")
                readBoundaryConditions(node, param_root, process_var);
            else if (node_name == "source_terms")
                readSourceTerms(node, param_root, process_var);
        }
        pvar = pvar.nextSibling();
    }
}

void XmlPrjInterface::readBoundaryConditions(
    QDomNode const& bc_root,
    QDomNode const& param_root,
    DataHolderLib::ProcessVariable const& pvar)
{
    QDomNode bc = bc_root.firstChild();
    while (bc != QDomNode())
    {
        std::unique_ptr<DataHolderLib::BoundaryCondition> cond(
            parseCondition<DataHolderLib::BoundaryCondition>(bc, param_root,
                                                             pvar));
        if (cond->getType() !=
            DataHolderLib::BoundaryCondition::ConditionType::NONE)
            _project.addBoundaryCondition(std::move(cond));

        bc = bc.nextSibling();
    }
}

void XmlPrjInterface::readSourceTerms(
    QDomNode const& st_root,
    QDomNode const& param_root,
    DataHolderLib::ProcessVariable const& pvar)
{
    QDomNode st = st_root.firstChild();
    while (st != QDomNode())
    {
        std::unique_ptr<DataHolderLib::SourceTerm> cond(
            parseCondition<DataHolderLib::SourceTerm>(st, param_root, pvar));
        if (cond->getType() != DataHolderLib::SourceTerm::ConditionType::NONE)
            _project.addSourceTerm(std::move(cond));
        st = st.nextSibling();
    }
}

template <typename T>
T* XmlPrjInterface::parseCondition(
    QDomNode const& node,
    QDomNode const& param_root,
    DataHolderLib::ProcessVariable const& pvar) const
{
    DataHolderLib::BaseObjType base_obj_type =
        DataHolderLib::BaseObjType::GEOMETRY;
    typename T::ConditionType type = T::ConditionType::NONE;
    std::string base_obj_name(""), obj_name(""), param_name("");
    QDomNode param_node = QDomNode();
    QDomNodeList nodeList = node.childNodes();
    for (int i = 0; i < nodeList.count(); i++)
    {
        QString const node_name = nodeList.at(i).nodeName();
        QString const content = nodeList.at(i).toElement().text().trimmed();
        if (node_name == "geometrical_set" &&
            base_obj_type != DataHolderLib::BaseObjType::MESH)
            base_obj_name = content.toStdString();
        else if (node_name == "geometry" &&
                 base_obj_type != DataHolderLib::BaseObjType::MESH)
            obj_name = content.toStdString();
        else if (node_name == "type")
            type = T::convertStringToType(content.toStdString());
        else if (node_name == "mesh")
        {
            base_obj_type = DataHolderLib::BaseObjType::MESH;
            base_obj_name = content.toStdString();
            obj_name.clear();
        }
        else if (node_name == "field_name")
            param_name = content.toStdString();
        else if (node_name == "parameter")
        {
            QDomNode val = findParam(param_root, content);
            if (val == QDomNode())
                continue;
            param_name = content.toStdString();
            param_node = nodeList.at(i);
        }
    }

    if (!param_name.empty())
    {
        T* cond = new T(pvar, param_name, type);
        if (base_obj_type == DataHolderLib::BaseObjType::MESH)
            cond->setMesh(base_obj_name);
        else
            cond->setGeoObject(base_obj_name, obj_name);

        return cond;
    }
    return new T({"", 0, 0}, "", T::ConditionType::NONE);
}

int XmlPrjInterface::writeToFile(const std::string& filename)
{
    _filename = filename;
    return BaseLib::IO::Writer::writeToFile(filename);
}

bool XmlPrjInterface::write()
{
    GeoLib::GEOObjects& geo_objects = _project.getGEOObjects();
    QFileInfo fi(QString::fromStdString(_filename));
    std::string path((fi.absolutePath()).toStdString() + "/");

    _out << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";  // xml
                                                                  // definition
    _out << "<?xml-stylesheet type=\"text/xsl\" "
            "href=\"OpenGeoSysProject.xsl\"?>\n\n";  // stylefile definition

    QDomDocument doc("OGS-PROJECT-DOM");
    QDomElement root = doc.createElement("OpenGeoSysProject");
    root.setAttribute("xmlns:ogs", "http://www.opengeosys.org");
    root.setAttribute("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    root.setAttribute(
        "xsi:noNamespaceSchemaLocation",
        "https://www.opengeosys.org/images/xsd/OpenGeoSysProject.xsd");

    doc.appendChild(root);

    // meshes
    auto const& mesh_vector = _project.getMeshObjects();
    for (auto const& mesh : mesh_vector)
    {
        // write mesh file
        MeshLib::IO::writeMeshToFile(*mesh, path + mesh->getName() + ".vtu");

        // write entry in project file
        QDomElement mesh_tag = doc.createElement("mesh");
        root.appendChild(mesh_tag);
        QDomText filename_text =
            doc.createTextNode(QString::fromStdString(mesh->getName()));
        mesh_tag.appendChild(filename_text);
    }

    // geometries
    std::vector<std::string> geo_names;
    geo_objects.getGeometryNames(geo_names);
    for (std::string const name : geo_names)
    {
        // write gml file
        GeoLib::IO::XmlGmlInterface gml(geo_objects);
        gml.setNameForExport(name);
        if (gml.writeToFile(std::string(path + name + ".gml")))
        {
            // write entry in project file
            QDomElement geo_tag = doc.createElement("geometry");
            root.appendChild(geo_tag);
            QDomText filename_text =
                doc.createTextNode(QString::fromStdString(name + ".gml"));
            geo_tag.appendChild(filename_text);
        }
        else
            ERR("XmlGmlInterface::writeFile(): Error writing gml-file \"%s\".",
                name.c_str());
    }

    // stations
    std::vector<std::string> stn_names;
    geo_objects.getStationVectorNames(stn_names);
    for (std::string const name : stn_names)
    {
        // write station file
        GeoLib::IO::XmlStnInterface stn(geo_objects);
        stn.setNameForExport(name);

        if (stn.writeToFile(path + name + ".stn"))
        {
            // write entry in project file
            QDomElement stn_tag = doc.createElement("stations");
            root.appendChild(stn_tag);
            QDomText filename_text =
                doc.createTextNode(QString::fromStdString(name + ".stn"));
            stn_tag.appendChild(filename_text);
        }
        else
            ERR("XmlStnInterface::writeFile(): Error writing stn-file \"%s\".",
                name.c_str());
    }

    if (_project.getBoundaryConditions().size() > 0 ||
        _project.getSourceTerms().size() > 0)
    {
        // parameters
        writeProcessVariables(doc, root);
    }

    std::string xml = doc.toString().toStdString();
    _out << xml;
    return true;
}

void addTextNode(QDomDocument& doc,
                 QDomElement& parent,
                 QString const& node_name,
                 QString const& content)
{
    QDomElement tag = doc.createElement(node_name);
    parent.appendChild(tag);
    QDomText order_text = doc.createTextNode(content);
    tag.appendChild(order_text);
}

bool PVarExists(std::string const& name,
                std::vector<DataHolderLib::ProcessVariable> const& p_vars)
{
    return std::any_of(p_vars.begin(), p_vars.end(),
                       [&](auto const& p_var) { return p_var.name == name; });
}

std::vector<DataHolderLib::ProcessVariable>
XmlPrjInterface::getPrimaryVariableVec() const
{
    std::vector<std::unique_ptr<DataHolderLib::BoundaryCondition>> const&
        boundary_conditions = _project.getBoundaryConditions();
    std::vector<std::unique_ptr<DataHolderLib::SourceTerm>> const&
        source_terms = _project.getSourceTerms();

    std::vector<DataHolderLib::ProcessVariable> p_vars;
    for (auto& bc : boundary_conditions)
    {
        DataHolderLib::ProcessVariable const& pvar(bc->getProcessVar());
        if (!PVarExists(pvar.name, p_vars))
            p_vars.push_back(pvar);
    }

    for (auto& st : source_terms)
    {
        DataHolderLib::ProcessVariable const& pvar(st->getProcessVar());
        if (!PVarExists(pvar.name, p_vars))
            p_vars.push_back(pvar);
    }
    return p_vars;
}

template <typename T>
void XmlPrjInterface::writeCondition(
    QDomDocument& doc,
    QDomElement& tag,
    DataHolderLib::FemCondition const& cond) const
{
    if (cond.getBaseObjType() == DataHolderLib::BaseObjType::GEOMETRY)
    {
        addTextNode(doc,
                    tag,
                    "geometrical_set",
                    QString::fromStdString(cond.getBaseObjName()));
        addTextNode(doc, tag, "geometry",
                    QString::fromStdString(cond.getObjName()));
        addTextNode(doc,
                    tag,
                    "type",
                    QString::fromStdString(T::convertTypeToString(
                        static_cast<T const&>(cond).getType())));
        addTextNode(doc, tag, "parameter",
                    QString::fromStdString(cond.getParamName()));
    }
    else if (cond.getBaseObjType() == DataHolderLib::BaseObjType::MESH)
    {
        addTextNode(doc,
                    tag,
                    "type",
                    QString::fromStdString(T::convertTypeToString(
                        static_cast<T const&>(cond).getType())));
        addTextNode(doc, tag, "mesh",
                    QString::fromStdString(cond.getBaseObjName()));
        addTextNode(doc,
                    tag,
                    "field_name",
                    QString::fromStdString(cond.getParamName()));
    }
}

void XmlPrjInterface::writeBoundaryConditions(QDomDocument& doc,
                                              QDomElement& bc_list_tag,
                                              std::string const& name) const
{
    std::vector<std::unique_ptr<DataHolderLib::BoundaryCondition>> const&
        boundary_conditions = _project.getBoundaryConditions();
    for (auto& bc : boundary_conditions)
    {
        if (bc->getProcessVarName() != name)
            continue;
        QDomElement bc_tag = doc.createElement("boundary_condition");
        bc_list_tag.appendChild(bc_tag);
        writeCondition<DataHolderLib::BoundaryCondition>(doc, bc_tag, *bc);
    }
}

void XmlPrjInterface::writeSourceTerms(QDomDocument& doc,
                                       QDomElement& st_list_tag,
                                       std::string const& name) const
{
    std::vector<std::unique_ptr<DataHolderLib::SourceTerm>> const&
        source_terms = _project.getSourceTerms();
    for (auto& st : source_terms)
    {
        if (st->getProcessVarName() != name)
            continue;
        QDomElement st_tag = doc.createElement("source_term");
        st_list_tag.appendChild(st_tag);
        writeCondition<DataHolderLib::SourceTerm>(doc, st_tag, *st);
    }
}

void XmlPrjInterface::writeProcessVariables(QDomDocument& doc,
                                            QDomElement& root) const
{
    std::vector<DataHolderLib::ProcessVariable> const p_vars(
        getPrimaryVariableVec());

    QDomElement param_list_tag = doc.createElement("parameters");
    root.appendChild(param_list_tag);

    QDomElement pvar_list_tag = doc.createElement("process_variables");
    root.appendChild(pvar_list_tag);

    for (DataHolderLib::ProcessVariable const p_var : p_vars)
    {
        QDomElement pvar_tag = doc.createElement("process_variable");
        pvar_list_tag.appendChild(pvar_tag);
        addTextNode(doc, pvar_tag, "name", QString::fromStdString(p_var.name));
        addTextNode(doc, pvar_tag, "order", QString::number(p_var.order));
        addTextNode(doc, pvar_tag, "components",
                    QString::number(p_var.components));
        QDomElement bc_list_tag = doc.createElement("boundary_conditions");
        pvar_tag.appendChild(bc_list_tag);
        writeBoundaryConditions(doc, bc_list_tag, p_var.name);
        QDomElement st_list_tag = doc.createElement("source_terms");
        pvar_tag.appendChild(st_list_tag);
        writeSourceTerms(doc, st_list_tag, p_var.name);
    }
}
}  // namespace FileIO
