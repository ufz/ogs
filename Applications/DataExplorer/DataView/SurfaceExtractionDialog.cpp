// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "SurfaceExtractionDialog.h"

#include <QDoubleValidator>

SurfaceExtractionDialog::SurfaceExtractionDialog(QDialog* parent)
    : QDialog(parent)
{
    setupUi(this);
    this->xNormalEdit->setValidator(
        new QDoubleValidator(-1, 1, 3, xNormalEdit));
    this->yNormalEdit->setValidator(
        new QDoubleValidator(-1, 1, 3, yNormalEdit));
    this->zNormalEdit->setValidator(
        new QDoubleValidator(-1, 1, 3, zNormalEdit));
}

void SurfaceExtractionDialog::accept()
{
    _dir = Eigen::Vector3d({xNormalEdit->text().toDouble(),
                            yNormalEdit->text().toDouble(),
                            zNormalEdit->text().toDouble()});
    _tolerance = degreesSpinBox->text().toInt();

    this->done(QDialog::Accepted);
}
