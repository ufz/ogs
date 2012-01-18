/**
 * \file ConditionWriterDialog.cpp
 * 2012/01/11 KR Initial implementation
 */

#include "ConditionWriterDialog.h"
#include "FEMCondition.h"
#include "OGSError.h"

#include <QFileDialog>
#include <QFileInfo>
#include <QSettings>

ConditionWriterDialog::ConditionWriterDialog(const GEOLIB::GEOObjects *geo_objects, QDialog* parent)
	: QDialog(parent)
{
	setupUi(this);

	std::vector<std::string> geo_names;
	geo_objects->getGeometryNames(geo_names);

	for (size_t i=0; i<geo_names.size(); i++)
		this->geoBox->addItem(QString::fromStdString(geo_names[i]));

}

ConditionWriterDialog::~ConditionWriterDialog()
{
}

void ConditionWriterDialog::on_fileNameButton_pressed()
{
	QString filetypes("");
	int geo_idx = this->geoBox->currentIndex();
	int cnd_idx = this->condTypeBox->currentIndex();
	if ((geo_idx == 0) || (cnd_idx == 0)) filetypes = "OpenGeoSys FEM Condition file (*.cnd)";
	else if ((geo_idx != 0) && (cnd_idx == 1)) filetypes = "OpenGeoSys FEM Condition file (*.cnd);;GeoSys Boundary Condition (*.bc)";
	else if ((geo_idx != 0) && (cnd_idx == 2)) filetypes = "OpenGeoSys FEM Condition file (*.cnd);;GeoSys Initial Condition (*.ic)";
	else if ((geo_idx != 0) && (cnd_idx == 3)) filetypes = "OpenGeoSys FEM Condition file (*.cnd);;GeoSys Source Condition (*.st)";

	QSettings settings("UFZ", "OpenGeoSys-5");
	QString fileName = QFileDialog::getSaveFileName(this, "Select path",
					settings.value("lastOpenedFileDirectory").toString(), filetypes);

	if (!fileName.isEmpty())
	{
		this->fileNameEdit->setText(fileName);

		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedFileDirectory", dir.absolutePath());
	}
}


void ConditionWriterDialog::accept()
{
	const QString file_name = this->fileNameEdit->text();
	QFileInfo fi(file_name);

	if ((fi.suffix().toLower().compare("cnd") != 0) &&
		((this->geoBox->currentIndex()==0) || (this->condTypeBox->currentIndex()==0)))
	{
		OGSError::box("Multiple geometries or multiple types of conditions\ncan only be saved to *.cnd files.","Inconsistent selection of parameters");
	}
	else
	{
		QString geo_name = this->geoBox->currentText();
		if (this->geoBox->currentIndex() == 0) geo_name = "";

		FEMCondition::CondType cond_type(FEMCondition::UNSPECIFIED);;
		switch (this->condTypeBox->currentIndex())
		{
			case 0: 
				cond_type = FEMCondition::UNSPECIFIED; break;
			case 1: 
				cond_type = FEMCondition::BOUNDARY_CONDITION; break;
			case 2: 
				cond_type = FEMCondition::INITIAL_CONDITION; break;
			case 3: 
				cond_type = FEMCondition::SOURCE_TERM; break;
			default:
				std::cout << "Error in ConditionWriterDialog..." << std::endl;
		}

		emit saveFEMConditionsRequested(geo_name, cond_type, file_name);

		this->done(QDialog::Accepted);
	}
}

void ConditionWriterDialog::reject()
{
	this->done(QDialog::Rejected);
}

