/**
 * \file OGSFileConverter.cpp
 * 2012/04/04 KR Initial implementation
 */
#define BOOST_FILESYSTEM_VERSION 3

#include "OGSFileConverter.h"
#include "FileListDialog.h"
#include "FileTools.h"
#include "OGSError.h"

#include <QFileInfo>

// conversion includes
#include "Applications/ApplicationsLib/ProjectData.h"
#include "StringTools.h"

// geometry
#include "GEOObjects.h"
#include "Legacy/OGSIOVer4.h"
#include "XmlIO/Qt/XmlGmlInterface.h"

// mesh
#include "Legacy/MeshIO.h"
#include "XmlIO/Boost/BoostVtuInterface.h"

OGSFileConverter::OGSFileConverter(QWidget* parent)
	: QDialog(parent)
{
	setupUi(this);
}

OGSFileConverter::~OGSFileConverter()
{
}

void OGSFileConverter::convertGML2GLI(const QStringList &input, const QString &output)
{
	if (input.empty())
		return;

	ProjectData project;
	GeoLib::GEOObjects* geo_objects = project.getGEOObjects();
	FileIO::XmlGmlInterface xml(*geo_objects);

	for (QStringList::const_iterator it=input.begin(); it!=input.end(); ++it)
	{
		const QFileInfo fi(*it);
		const std::string file_name = fi.baseName().toStdString();
		const std::string output_str = QString(output + "/" + fi.completeBaseName() + ".gli").toStdString();

		if (!fileExists(output_str))
		{
			xml.readFile(*it);
			std::vector<std::string> geo_names;
			geo_objects->getGeometryNames(geo_names);
			FileIO::Legacy::writeGLIFileV4(output_str, geo_names[0], *geo_objects);
			geo_objects->removeSurfaceVec(geo_names[0]);
			geo_objects->removePolylineVec(geo_names[0]);
			geo_objects->removePointVec(geo_names[0]);
		}
	}

	//FileIO::writeAllDataToGLIFileV4(output.toStdString(), *geo_objects);
	OGSError::box("File conversion finished");
}

void OGSFileConverter::convertGLI2GML(const QStringList &input, const QString &output)
{
	if (input.empty())
		return;

	ProjectData project;
	GeoLib::GEOObjects* geo_objects = project.getGEOObjects();
	FileIO::XmlGmlInterface xml(*geo_objects);

	for (QStringList::const_iterator it=input.begin(); it!=input.end(); ++it)
	{
		const QFileInfo fi(*it);
		const std::string output_str = QString(output + "/" + fi.completeBaseName() + ".gml").toStdString();
		const std::string geo_name = BaseLib::extractBaseName(it->toStdString());

		if (!fileExists(output_str))
		{
			std::string unique_name;
			std::vector<std::string> errors;

			FileIO::Legacy::readGLIFileV4(it->toStdString(), geo_objects, unique_name, errors);
			if (errors.empty() || (errors.size()==1 && errors[0].compare("[readSurface] polyline for surface not found!")==0))
			{
				xml.setNameForExport(geo_name);
				xml.writeToFile(output_str);
				geo_objects->removeSurfaceVec(geo_name);
				geo_objects->removePolylineVec(geo_name);
				geo_objects->removePointVec(geo_name);
			}
			else
				for (size_t k(0); k<errors.size(); k++)
					OGSError::box(QString::fromStdString(errors[k]));
		}
	}
	
	OGSError::box("File conversion finished");
}

void OGSFileConverter::convertVTU2MSH(const QStringList &input, const QString &output)
{
	if (input.empty())
		return;

	for (QStringList::const_iterator it=input.begin(); it!=input.end(); ++it)
	{
		const QFileInfo fi(*it);
		const std::string msh_name = fi.fileName().toStdString();
		const std::string output_str = QString(output + "/" + fi.completeBaseName() + ".msh").toStdString();

		if (!fileExists(output_str))
		{
			FileIO::BoostVtuInterface vtu;
			MeshLib::Mesh const*const mesh (vtu.readVTUFile(it->toStdString().c_str()));
			FileIO::Legacy::MeshIO meshIO;
			meshIO.setMesh(mesh);
			meshIO.writeToFile(output_str.c_str());
			delete mesh;
		}
	}

	OGSError::box("File conversion finished");
}

void OGSFileConverter::convertMSH2VTU(const QStringList &input, const QString &output)
{
	if (input.empty())
		return;

	for (QStringList::const_iterator it=input.begin(); it!=input.end(); ++it)
	{
		const QFileInfo fi(*it);
		const std::string output_str = QString(output + "/" + fi.completeBaseName() + ".vtu").toStdString();
		const std::string msh_name = BaseLib::extractBaseName(it->toStdString());

		if (!fileExists(output_str))
		{
			FileIO::Legacy::MeshIO meshIO;
			MeshLib::Mesh const*const mesh (meshIO.loadMeshFromFile(it->toStdString()));
			FileIO::BoostVtuInterface vtu;
			vtu.setMesh(mesh);
			vtu.writeToFile(output_str);
			delete mesh;
		}
	}
	
	OGSError::box("File conversion finished");
}

void OGSFileConverter::on_gml2gliButton_pressed()
{
	FileListDialog dlg(FileListDialog::GML, FileListDialog::GLI);
	if (dlg.exec())
		convertGML2GLI(dlg.getInputFileList(), dlg.getOutputDir());
}

void OGSFileConverter::on_gli2gmlButton_pressed()
{
	FileListDialog dlg(FileListDialog::GLI, FileListDialog::GML);
	if (dlg.exec())
		convertGLI2GML(dlg.getInputFileList(), dlg.getOutputDir());
}

void OGSFileConverter::on_vtu2mshButton_pressed()
{
	FileListDialog dlg(FileListDialog::VTU, FileListDialog::MSH);
	if (dlg.exec())
		convertVTU2MSH(dlg.getInputFileList(), dlg.getOutputDir());
}

void OGSFileConverter::on_msh2vtuButton_pressed()
{
	FileListDialog dlg(FileListDialog::MSH, FileListDialog::VTU);
	if (dlg.exec())
		convertMSH2VTU(dlg.getInputFileList(), dlg.getOutputDir());
}

void OGSFileConverter::on_closeDialogButton_pressed()
{
	this->close();
}

bool OGSFileConverter::fileExists(const std::string &file_name) const
{
	std::ifstream file(file_name.c_str());
	if (file)
	{
		QString name = QString::fromStdString(BaseLib::extractBaseName(file_name));
		return !OGSError::question("The file \'" + name + "\' already exists.\n Do you want to overwrite it?", "Warning");
	}
	return false;
}

