/**
 * \file CondFromRasterDialog.cpp
 * 2012/01/04 KR Initial implementation
 */

#include "CondFromRasterDialog.h"

#include <QFileDialog>
#include <QSettings>

#include "DirectConditionGenerator.h"

CondFromRasterDialog::CondFromRasterDialog(const ProjectData *project, QDialog* parent)
	: QDialog(parent), _project(project)
{
	setupUi(this);

	const std::map<std::string, MeshLib::CFEMesh*> msh_map = project->getMeshObjects();

	for (std::map<std::string, MeshLib::CFEMesh*>::const_iterator it = msh_map.begin();
			                                                      it != msh_map.end(); ++it)
		 this->meshBox->addItem(QString::fromStdString(it->first));

	this->directButton->setChecked(true);
}

CondFromRasterDialog::~CondFromRasterDialog()
{
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
	const MeshLib::CFEMesh* mesh = _project->getMesh(mesh_name);

	if (this->directButton->isChecked())
	{
		DirectConditionGenerator dcg;
		dcg.fromRasterToSurfaceNodes(*mesh, raster_name);
		dcg.writeToFile(raster_name + ".txt");
	}
	else
	{
		// send integrate signal
	}
	this->done(QDialog::Accepted);
}

void CondFromRasterDialog::reject()
{
	this->done(QDialog::Rejected);
}

