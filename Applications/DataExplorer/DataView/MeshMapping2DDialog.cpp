/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshMapping2DDialog.h"
#include "OGSError.h"
#include "StrictDoubleValidator.h"

#include <QSettings>
#include <QFileDialog>

MeshMapping2DDialog::MeshMapping2DDialog(QDialog* parent)
: QDialog(parent)
{
    setupUi(this);

    StrictDoubleValidator* no_data_validator = new StrictDoubleValidator(this);
    this->noDataValueEdit->setValidator (no_data_validator);
    StrictDoubleValidator* static_value_validator = new StrictDoubleValidator(this);
    this->staticValueEdit->setValidator (static_value_validator);
}

void MeshMapping2DDialog::on_rasterValueButton_toggled(bool isChecked)
{
    this->rasterPathEdit->setEnabled(isChecked);
    this->noDataValueEdit->setEnabled(isChecked);
    this->rasterSelectButton->setEnabled(isChecked);
    this->staticValueEdit->setEnabled(!isChecked);
}

void MeshMapping2DDialog::on_rasterSelectButton_pressed()
{
    QSettings settings;
    QString filename = QFileDialog::getOpenFileName(
        this, "Select raster file to open",
        settings.value("lastOpenedRasterFileDirectory").toString(),
        "ASCII raster files (*.asc);;All files (* *.*)");
    this->rasterPathEdit->setText(filename);
    QFileInfo fi(filename);
    settings.setValue("lastOpenedRasterFileDirectory", fi.absolutePath());
}

void MeshMapping2DDialog::accept()
{
    if (this->rasterValueButton->isChecked() && this->rasterPathEdit->text().isEmpty())
    {
        OGSError::box("Please specify path to raster file.");
        return;
    }
    if (this->rasterValueButton->isChecked() && this->noDataValueEdit->text().isEmpty())
    {
        OGSError::box("Please specify No Data value.");
        return;
    }
    if (this->staticValueButton->isChecked() && this->staticValueEdit->text().isEmpty())
    {
        OGSError::box("Please specify value for mapping.");
        return;
    }
    if (this->newNameEdit->text().isEmpty())
    {
        OGSError::box("Please specify a name for the resulting mesh.");
        return;
    }

    this->done(QDialog::Accepted);
}

