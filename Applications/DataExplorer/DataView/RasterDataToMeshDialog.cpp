/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RasterDataToMeshDialog.h"

#include <QFileDialog>
#include <QSettings>

#include "Base/OGSError.h"
#include "Base/StrictDoubleValidator.h"

RasterDataToMeshDialog::RasterDataToMeshDialog(std::string const& mesh_name,
                                               QDialog* parent)
    : QDialog(parent)
{
    setupUi(this);

    this->meshNameEdit->setText(QString::fromStdString(mesh_name) + "_Raster");
    auto* no_data_validator = new StrictDoubleValidator(this);
    this->noDataValueEdit->setValidator(no_data_validator);
}

void RasterDataToMeshDialog::on_rasterSelectButton_pressed()
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

void RasterDataToMeshDialog::accept()
{
    if (this->rasterPathEdit->text().isEmpty())
    {
        OGSError::box("Please specify path to raster file.");
        return;
    }
    else if (this->arrayNameEdit->text().isEmpty())
    {
        OGSError::box("Please specify a name for the new array.");
        return;
    }
    else if (this->noDataValueEdit->text().isEmpty())
    {
        OGSError::box("Please specify No Data value.");
        return;
    }
    else if (this->meshNameEdit->text().isEmpty())
    {
        OGSError::box("Please specify name of new mesh.");
        return;
    }

    this->done(QDialog::Accepted);
}
