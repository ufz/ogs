/**
 * \file CondFromRasterDialog.cpp
 * 2012/01/04 KR Initial implementation
 */

#include "CondFromRasterDialog.h"

#include <QFileDialog>
#include <QSettings>

#include "DirectConditionGenerator.h"
#include "OGSError.h"
#include "StrictDoubleValidator.h"

CondFromRasterDialog::CondFromRasterDialog(const ProjectData *project, QDialog* parent)
	: QDialog(parent), _project(project)
{
	setupUi(this);

	this->scalingEdit->setEnabled(false);
	_scale_validator = new StrictDoubleValidator(-1e+10, 1e+10, 5);
	this->scalingEdit->setText("1.0");
	this->scalingEdit->setValidator (_scale_validator);

	const std::map<std::string, MeshLib::CFEMesh*> msh_map = project->getMeshObjects();

	for (std::map<std::string, MeshLib::CFEMesh*>::const_iterator it = msh_map.begin();
			                                                      it != msh_map.end(); ++it)
		 this->meshBox->addItem(QString::fromStdString(it->first));

	this->directButton->setChecked(true);
}

CondFromRasterDialog::~CondFromRasterDialog()
{
	delete _scale_validator;
}

void CondFromRasterDialog::on_selectButton_pressed()
{
	QSettings settings("UFZ", "OpenGeoSys-5");
#ifdef libgeotiff_FOUND
	QString geotiffExtension(" *.tif");
#else
	QString geotiffExtension("");
#endif
	QString fileName = QFileDialog::getOpenFileName(this, "Select raster file",
					settings.value("lastOpenedFileDirectory").toString(), QString(
									"Raster files (*.asc *.bmp *.jpg *.png%1);;") .arg(geotiffExtension));

	if (!fileName.isEmpty())
	{
		this->rasterEdit->setText(fileName);

		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
}


void CondFromRasterDialog::accept()
{
	std::string mesh_name (this->meshBox->currentText().toStdString());
	std::string raster_name (this->rasterEdit->text().toStdString());
	double scaling_factor = strtod(this->scalingEdit->text().toStdString().c_str(), 0);

	if (mesh_name.empty())
	{
		OGSError::box("No mesh selected.");
		return;
	}
	if (raster_name.empty())
	{
		OGSError::box("No raster selected.");
		return;
	}
	
	const MeshLib::CFEMesh* mesh = _project->getMesh(mesh_name);

	if (this->directButton->isChecked())
	{
		DirectConditionGenerator dcg;
		dcg.directToSurfaceNodes(*mesh, raster_name);
		dcg.writeToFile(raster_name + ".txt");
	}
	else
	{
		if (scaling_factor == 0)
		{
			OGSError::box("No valid scaling factor given.");
			return;
		}
		MeshLib::CFEMesh* new_mesh = const_cast<MeshLib::CFEMesh*>(mesh);
		DirectConditionGenerator dcg;
		dcg.directWithSurfaceIntegration(*new_mesh, raster_name);
		dcg.writeToFile(raster_name + ".txt");
	}
	this->done(QDialog::Accepted);
}

void CondFromRasterDialog::reject()
{
	this->done(QDialog::Rejected);
}

void CondFromRasterDialog::on_integrateButton_toggled(bool isSelected)
{
	this->scalingEdit->setEnabled(isSelected);
}