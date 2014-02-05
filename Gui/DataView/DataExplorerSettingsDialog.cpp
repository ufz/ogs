/**
 * \file   DataExplorerSettingsDialog.cpp
 * \author Karsten Rink
 * \date   2014-02-05
 * \brief  Implementation of the DataExplorerSettingsDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DataExplorerSettingsDialog.h"

#include <QFileDialog>
#include <QSettings>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "OGSError.h"


DataExplorerSettingsDialog::DataExplorerSettingsDialog(QDialog* parent) : QDialog(parent)
{
	setupUi(this);

	QSettings settings("UFZ", "OpenGeoSys");
	this->fcPathEdit->setText(settings.value("DataExplorerConverterPath").toString());
	this->gmshPathEdit->setText(settings.value("DataExplorerGMSHPath").toString());
}

DataExplorerSettingsDialog::~DataExplorerSettingsDialog()
{
}

void DataExplorerSettingsDialog::on_fcPathButton_clicked()
{
	QSettings settings("UFZ", "OpenGeoSys");
	QString file_name = QFileDialog::getOpenFileName(this, "Select path for OGS File Converter...",
						                                   settings.value("DataExplorerConverterPath").toString(),
						                                   "OGSFileConverter.exe;;");

	if (!file_name.isEmpty())
		this->fcPathEdit->setText(file_name);
}

void DataExplorerSettingsDialog::on_gmshPathButton_clicked()
{
	QSettings settings("UFZ", "OpenGeoSys");
	QString file_name = QFileDialog::getOpenFileName(this, "Select path for GMSH...",
						                                   settings.value("DataExplorerGMSHPath").toString(),
						                                   "gmsh*;;");

	if (!file_name.isEmpty())
		this->gmshPathEdit->setText(file_name);
}

void DataExplorerSettingsDialog::accept()
{
	QSettings settings("UFZ", "OpenGeoSys");
	settings.setValue("DataExplorerConverterPath", this->fcPathEdit->text());
	settings.setValue("DataExplorerGMSHPath", this->gmshPathEdit->text());
	this->done(QDialog::Accepted);
}


