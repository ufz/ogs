/**
 * \file
 * \author Karsten Rink
 * \date   2012-01-11
 * \brief  Implementation of the ConditionWriterDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ConditionWriterDialog.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "FEMCondition.h"
#include "OGSError.h"
#include "LastSavedFileDirectory.h"

#include <QFileDialog>
#include <QFileInfo>
#include <QSettings>

ConditionWriterDialog::ConditionWriterDialog(const GeoLib::GEOObjects *geo_objects, QDialog* parent)
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
	QSettings settings;
	QString fileName = QFileDialog::getSaveFileName(this, "Select path",
					LastSavedFileDirectory::getDir(), "OpenGeoSys FEM Condition file (*.cnd)");

	if (!fileName.isEmpty())
		this->fileNameEdit->setText(fileName);
}

void ConditionWriterDialog::accept()
{
	const QString file_name = this->fileNameEdit->text();
	QFileInfo fi(file_name);

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
			ERR("ConditionWriterDialog::accept(): case %d not handled.", this->condTypeBox->currentIndex());
	}

	LastSavedFileDirectory::setDir(file_name);
	emit saveFEMConditionsRequested(geo_name, cond_type, file_name);

	this->done(QDialog::Accepted);
}

void ConditionWriterDialog::reject()
{
	this->done(QDialog::Rejected);
}

