/**
 * \file OGSFileConverter.cpp
 * 2012/04/04 KR Initial implementation
 */

#include "OGSFileConverter.h"
#include "FileListDialog.h"
#include "ConversionTools.h"
#include "OGSError.h"

#include <QFileInfo>

// conversion includes
#include "ProjectData.h"
#include "GEOObjects.h"
#include "OGSIOVer4.h"
#include "XmlIO/XmlCndInterface.h"
#include "XmlIO/XmlGmlInterface.h"
#include "StringTools.h"

// old condition objects
#include "BoundaryCondition.h"
#include "InitialCondition.h"
#include "SourceTerm.h"
#include "rf_bc_new.h"
#include "rf_ic_new.h"
#include "rf_st_new.h"

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
	ProjectData project;
	GeoLib::GEOObjects* geo_objects = new GeoLib::GEOObjects;
	project.setGEOObjects(geo_objects);

	FileFinder fileFinder = createFileFinder();
	std::string schemaName(fileFinder.getPath("OpenGeoSysGLI.xsd"));
	FileIO::XmlGmlInterface xml(&project, schemaName);

	for (QStringList::const_iterator it=input.begin(); it!=input.end(); ++it)
		xml.readFile(*it);

	FileIO::writeAllDataToGLIFileV4(output.toStdString(), *geo_objects);
	OGSError::box("File conversion finished");
}

void OGSFileConverter::convertGLI2GML(const QStringList &input, const QString &output)
{
	ProjectData project;
	GeoLib::GEOObjects* geo_objects = new GeoLib::GEOObjects;
	project.setGEOObjects(geo_objects);

	std::vector<std::string> merge_list;
	for (QStringList::const_iterator it=input.begin(); it!=input.end(); ++it)
	{
		std::string unique_name;
		std::vector<std::string> errors;
		if (! FileIO::readGLIFileV4(it->toStdString(), geo_objects, unique_name, errors))
		{
			for (size_t k(0); k<errors.size(); k++)
				OGSError::box(QString::fromStdString(errors[k]));
		}
		else
			merge_list.push_back(unique_name);
	}

	if (!merge_list.empty())
	{
		std::string merged_geo_name (merge_list[0]);
		if (merge_list.size()>1)
		{
			merged_geo_name = BaseLib::getFileNameFromPath(output.toStdString());
			geo_objects->mergeGeometries(merge_list, merged_geo_name);
		}
		FileFinder fileFinder = createFileFinder();
		std::string schemaName(fileFinder.getPath("OpenGeoSysGLI.xsd"));
		FileIO::XmlGmlInterface xml(&project, schemaName);
		xml.setNameForExport(merged_geo_name);
		xml.writeToFile(output.toStdString());
	}
	OGSError::box("File conversion finished");
}

void OGSFileConverter::convertCND2BC(const QStringList &input, const QString &output)
{
	ProjectData project;
	GeoLib::GEOObjects* geo_objects = new GeoLib::GEOObjects;
	project.setGEOObjects(geo_objects);

	// HACK for enabling conversion of files without loading the associated geometry
	std::vector<GeoLib::Point*> *fake_geo = new std::vector<GeoLib::Point*>;
	fake_geo->push_back(new GeoLib::Point(0,0,0));
	std::string fake_name("conversionTestRun#1");
	geo_objects->addPointVec(fake_geo, fake_name);

	FileFinder fileFinder = createFileFinder();
	std::string schemaName(fileFinder.getPath("OpenGeoSysCond.xsd"));
	FileIO::XmlCndInterface xml(&project, schemaName);

	std::vector<FEMCondition*> conditions;

	for (QStringList::const_iterator it=input.begin(); it!=input.end(); ++it)
		xml.readFile(conditions, *it);

	if (!conditions.empty())
	{
		project.addConditions(conditions);
		QFileInfo fi(output);
		FEMCondition::CondType type = FEMCondition::UNSPECIFIED;
		if (fi.suffix().compare("bc") == 0)      type = FEMCondition::BOUNDARY_CONDITION;
		else if (fi.suffix().compare("ic") == 0) type = FEMCondition::INITIAL_CONDITION;
		else if (fi.suffix().compare("st") == 0) type = FEMCondition::SOURCE_TERM;

		size_t count(0);
		size_t nConds(conditions.size());
		for (size_t i=0; i<nConds; i++)
		{
			if (conditions[i]->getCondType() == type)
			{
				if (type == FEMCondition::BOUNDARY_CONDITION)
					bc_list.push_back(new CBoundaryCondition(static_cast<BoundaryCondition*>(conditions[i])));
				else if (type == FEMCondition::INITIAL_CONDITION)
					ic_vector.push_back(new CInitialCondition(static_cast<InitialCondition*>(conditions[i])));
				else if (type == FEMCondition::SOURCE_TERM)
					st_vector.push_back(new CSourceTerm(static_cast<SourceTerm*>(conditions[i])));

				if (conditions[i]->getProcessDistributionType() == FiniteElement::DIRECT)
				{
					std::string count_str (QString::number(count++).toStdString());
					std::string direct_value_file = fi.absolutePath().toStdString() + "/direct_values" + count_str + ".txt";
					st_vector[st_vector.size()-1]->fname = direct_value_file;
					ConversionTools::writeDirectValues(*conditions[i], direct_value_file);
				}
			}
		}
		if (type == FEMCondition::BOUNDARY_CONDITION)
			BCWrite(output.toStdString());
		else if (type == FEMCondition::INITIAL_CONDITION)
			ICWrite(output.toStdString());
		else if (type == FEMCondition::SOURCE_TERM)
			STWrite(output.toStdString());
	}
	OGSError::box("File conversion finished");
}

void OGSFileConverter::convertBC2CND(const QStringList &input, const QString &output)
{
	ProjectData project;
	std::vector<FEMCondition*> conditions;
	for (QStringList::const_iterator it=input.begin(); it!=input.end(); ++it)
		ConversionTools::getFEMConditionsFromASCIIFile(*it, conditions);

	if (!conditions.empty())
	{
		project.addConditions(conditions);
		FileFinder fileFinder = createFileFinder();
		std::string schemaName(fileFinder.getPath("OpenGeoSysCond.xsd"));
		FileIO::XmlCndInterface xml(&project, schemaName);
		xml.writeToFile(output.toStdString());
	}
	OGSError::box("File conversion finished");
}

FileFinder OGSFileConverter::createFileFinder()
{
	FileFinder fileFinder;
	fileFinder.addDirectory(".");
	fileFinder.addDirectory(std::string(SOURCEPATH).append("/FileIO"));
	return fileFinder;
}

void OGSFileConverter::on_gml2gliButton_pressed()
{
	FileListDialog dlg(FileListDialog::GML, FileListDialog::GLI);
	connect(&dlg, SIGNAL(fileLists(const QStringList, const QString)),
	        this, SLOT(convertGML2GLI(const QStringList, const QString)));
	dlg.exec();
}

void OGSFileConverter::on_gli2gmlButton_pressed()
{
	FileListDialog dlg(FileListDialog::GLI, FileListDialog::GML);
	connect(&dlg, SIGNAL(fileLists(const QStringList, const QString)),
	        this, SLOT(convertGLI2GML(const QStringList, const QString)));
	dlg.exec();
}

void OGSFileConverter::on_bc2cndButton_pressed()
{
	FileListDialog dlg(FileListDialog::BC, FileListDialog::CND);
	connect(&dlg, SIGNAL(fileLists(const QStringList, const QString)),
	        this, SLOT(convertBC2CND(const QStringList, const QString)));
	dlg.exec();
}

void OGSFileConverter::on_cnd2bcButton_pressed()
{
	FileListDialog dlg(FileListDialog::CND, FileListDialog::BC);
	connect(&dlg, SIGNAL(fileLists(const QStringList, const QString)),
	        this, SLOT(convertCND2BC(const QStringList, const QString)));
	dlg.exec();
}

void OGSFileConverter::on_closeDialogButton_pressed()
{
	this->close();
}

