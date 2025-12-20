// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "DataExplorerSettingsDialog.h"

#include <QFileDialog>
#include <QSettings>

#include "Base/OGSError.h"

DataExplorerSettingsDialog::DataExplorerSettingsDialog(QDialog* parent)
    : QDialog(parent)
{
    setupUi(this);

    QSettings settings;
    this->gmshPathEdit->setText(
        settings.value("DataExplorerGmshPath").toString());
}

DataExplorerSettingsDialog::~DataExplorerSettingsDialog() = default;

void DataExplorerSettingsDialog::on_gmshPathButton_clicked()
{
    QSettings settings;
    QString file_name = QFileDialog::getOpenFileName(
        this, "Select path for GMSH...",
        settings.value("DataExplorerGmshPath").toString(), "*gmsh*");
    if (!file_name.isEmpty())
    {
        this->gmshPathEdit->setText(file_name);
    }
}

void DataExplorerSettingsDialog::accept()
{
    QSettings settings;
    settings.setValue("DataExplorerGmshPath", this->gmshPathEdit->text());
    this->done(QDialog::Accepted);
}
