/**
 * \file   OGSFileConverter.cpp
 * \author Karsten Rink
 * \date   2012-04-04
 * \brief  Implementation of OGSFileConverter class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#define BOOST_FILESYSTEM_VERSION 3

#include "OGSFileConverter.h"

#include <QFileInfo>

#include "Applications/DataExplorer/Base/OGSError.h"
#include "Applications/FileIO/Legacy/OGSIOVer4.h"
#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"
#include "FileListDialog.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/IO/XmlIO/Qt/XmlGmlInterface.h"
#include "MeshLib/IO/Legacy/MeshIO.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/Mesh.h"

OGSFileConverter::OGSFileConverter(QWidget* parent)
    : QDialog(parent)
{
    setupUi(this);
}

OGSFileConverter::~OGSFileConverter() = default;

void OGSFileConverter::convertGML2GLI(const QStringList &input, const QString &output) const
{
    if (input.empty())
        return;

    GeoLib::GEOObjects geo_objects;
    GeoLib::IO::XmlGmlInterface xml(geo_objects);

    for (const auto& input_string : input)
    {
        const QFileInfo fi(input_string);
        const std::string output_str = QString(output + "/" + fi.completeBaseName() + ".gli").toStdString();

        if (fileExists(output_str))
            continue;

        if (!xml.readFile(input_string))
        {
            OGSError::box("Error reading geometry " + fi.fileName());
            continue;
        }
        std::vector<std::string> geo_names;
        geo_objects.getGeometryNames(geo_names);
        FileIO::Legacy::writeGLIFileV4(output_str, geo_names[0], geo_objects);
        geo_objects.removeSurfaceVec(geo_names[0]);
        geo_objects.removePolylineVec(geo_names[0]);
        geo_objects.removePointVec(geo_names[0]);
    }
    OGSError::box("File conversion finished");
}

void OGSFileConverter::convertGLI2GML(const QStringList &input, const QString &output) const
{
    if (input.empty())
        return;

    GeoLib::GEOObjects geo_objects;
    GeoLib::IO::XmlGmlInterface xml(geo_objects);

    for (const auto& input_string : input)
    {
        const QFileInfo fi(input_string);
        const std::string output_str = QString(output + "/" + fi.completeBaseName() + ".gml").toStdString();

        if (fileExists(output_str))
            continue;

        std::string unique_name;
        std::vector<std::string> errors;

        FileIO::Legacy::readGLIFileV4(input_string.toStdString(),
                                          geo_objects, unique_name, errors);
        if (errors.empty() ||
            (errors.size() == 1 &&
             errors[0] == "[readSurface] polyline for surface not found!"))
        {
            std::string const geo_name =
                BaseLib::extractBaseName(input_string.toStdString());
            xml.setNameForExport(geo_name);
            xml.writeToFile(output_str);
            geo_objects.removeSurfaceVec(geo_name);
            geo_objects.removePolylineVec(geo_name);
            geo_objects.removePointVec(geo_name);
        }
        else
            for (auto& error : errors)
                OGSError::box(QString::fromStdString(error));
    }
    OGSError::box("File conversion finished");
}

void OGSFileConverter::convertVTU2MSH(const QStringList &input, const QString &output) const
{
    if (input.empty())
        return;

    for (const auto& input_string : input)
    {
        const QFileInfo fi(input_string);
        const std::string output_str = QString(output + "/" + fi.completeBaseName() + ".msh").toStdString();

        if (fileExists(output_str))
            continue;

        MeshLib::Mesh const* const mesh(MeshLib::IO::VtuInterface::readVTUFile(
            input_string.toStdString().c_str()));
        if (mesh == nullptr)
        {
            OGSError::box("Error reading mesh " + fi.fileName());
            continue;
        }
        MeshLib::IO::Legacy::MeshIO meshIO;
        meshIO.setMesh(mesh);
        meshIO.writeToFile(output_str.c_str());
        delete mesh;
    }
    OGSError::box("File conversion finished");
}

void OGSFileConverter::convertMSH2VTU(const QStringList &input, const QString &output) const
{
    if (input.empty())
        return;

    for (const auto& input_string : input)
    {
        const QFileInfo fi(input_string);
        const std::string output_str = QString(output + "/" + fi.completeBaseName() + ".vtu").toStdString();

        if (fileExists(output_str))
            continue;

        MeshLib::IO::Legacy::MeshIO meshIO;
        MeshLib::Mesh const* const mesh(
            meshIO.loadMeshFromFile(input_string.toStdString()));
        if (mesh == nullptr)
        {
            OGSError::box("Error reading mesh " + fi.fileName());
            continue;
        }
        MeshLib::IO::VtuInterface vtu(mesh);
        vtu.writeToFile(output_str);
        delete mesh;
    }
    OGSError::box("File conversion finished");
}

void OGSFileConverter::on_gml2gliButton_pressed() const
{
    FileListDialog dlg(FileType::GML, FileType::GLI);
    if (dlg.exec())
        convertGML2GLI(dlg.getInputFileList(), dlg.getOutputDir());
}

void OGSFileConverter::on_gli2gmlButton_pressed() const
{
    FileListDialog dlg(FileType::GLI, FileType::GML);
    if (dlg.exec())
        convertGLI2GML(dlg.getInputFileList(), dlg.getOutputDir());
}

void OGSFileConverter::on_vtu2mshButton_pressed() const
{
    FileListDialog dlg(FileType::VTU, FileType::MSH);
    if (dlg.exec())
        convertVTU2MSH(dlg.getInputFileList(), dlg.getOutputDir());
}

void OGSFileConverter::on_msh2vtuButton_pressed() const
{
    FileListDialog dlg(FileType::MSH, FileType::VTU);
    if (dlg.exec())
        convertMSH2VTU(dlg.getInputFileList(), dlg.getOutputDir());
}

void OGSFileConverter::on_closeDialogButton_pressed()
{
    this->close();
}

bool OGSFileConverter::fileExists(const std::string &file_name) const
{
    std::ifstream const file(file_name.c_str());
    if (file)
    {
        QString const name = QString::fromStdString(BaseLib::extractBaseName(file_name));
        return !OGSError::question("The file \'" + name + "\' already exists.\n Do you want to overwrite it?", "Warning");
    }
    return false;
}
