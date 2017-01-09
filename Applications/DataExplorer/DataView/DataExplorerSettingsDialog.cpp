/**
 * \file   DataExplorerSettingsDialog.cpp
 * \author Karsten Rink
 * \date   2014-02-05
 * \brief  Implementation of the DataExplorerSettingsDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DataExplorerSettingsDialog.h"

#include <QFileDialog>
#include <QSettings>

#include <logog/include/logog.hpp>

#include "OGSError.h"


DataExplorerSettingsDialog::DataExplorerSettingsDialog(QDialog* parent) : QDialog(parent)
{
    setupUi(this);

    QSettings settings;
    this->gmshPathEdit->setText(settings.value("DataExplorerGmshPath").toString());
}

DataExplorerSettingsDialog::~DataExplorerSettingsDialog()
{
}

void DataExplorerSettingsDialog::on_gmshPathButton_clicked()
{
    QSettings settings;
    QString file_name = QFileDialog::getOpenFileName(this, "Select path for GMSH...",
                                                           settings.value("DataExplorerGmshPath").toString(),
                                                           "*gmsh*");
    if (!file_name.isEmpty())
        this->gmshPathEdit->setText(file_name);
}

void DataExplorerSettingsDialog::accept()
{
    QSettings settings;
    settings.setValue("DataExplorerGmshPath", this->gmshPathEdit->text());
    this->done(QDialog::Accepted);
}


