/**
 * \file
 * \author Karsten Rink
 * \date   2011-03-16
 * \brief  Implementation of the MshQualitySelectionDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshQualitySelectionDialog.h"
#include "OGSError.h"

#include <QFileDialog>
#include <QFileInfo>
#include <QSettings>

/// Constructor
MeshQualitySelectionDialog::MeshQualitySelectionDialog(QDialog* parent)
: QDialog(parent), _metric (MeshQualityType::EDGERATIO), _histogram_path("")
{
	setupUi(this);
	this->choiceEdges->toggle();
}

MeshQualitySelectionDialog::~MeshQualitySelectionDialog()
{
}

void MeshQualitySelectionDialog::on_histogramCheckBox_toggled(bool is_checked) const
{
	histogramPathEdit->setEnabled(is_checked); 
	histogramPathButton->setEnabled(is_checked);
}

void MeshQualitySelectionDialog::on_histogramPathButton_pressed()
{
	QSettings settings;
	QFileInfo fi(settings.value("lastOpenedFileDirectory").toString());
	QString file_name = QFileDialog::getSaveFileName(this, "Save histogram as",
		fi.baseName(),
		"Text files (*.txt);;All files (* *.*)");
	this->histogramPathEdit->setText(file_name);
}

/// Instructions if the OK-Button has been pressed.
void MeshQualitySelectionDialog::accept()
{
	if (this->choiceEdges->isChecked())
		_metric = MeshQualityType::EDGERATIO;
	else if (this->choiceArea->isChecked())
		_metric = MeshQualityType::ELEMENTSIZE;
	else if (this->choiceVolume->isChecked())
		_metric = MeshQualityType::SIZEDIFFERENCE;
	else if (this->choiceAngles->isChecked())
		_metric = MeshQualityType::EQUIANGLESKEW;
	else if (this->choiceRadius->isChecked())
		_metric = MeshQualityType::RADIUSEDGERATIO;
	else
		_metric = MeshQualityType::INVALID;

	if (this->histogramCheckBox->isChecked())
	{
		_histogram_path = this->histogramPathEdit->text().toStdString();
		if (_histogram_path.empty())
		{
			OGSError::box("No path for histogram file specified.");
			return;
		}
	}

	this->done(QDialog::Accepted);
}

/// Instructions if the Cancel-Button has been pressed.
void MeshQualitySelectionDialog::reject()
{
	this->done(QDialog::Rejected);
}
