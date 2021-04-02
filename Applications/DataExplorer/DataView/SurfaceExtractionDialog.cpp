/**
 * \file
 * \author Karsten Rink
 * \date   2015-01-29
 * \brief  Implementation of the SaveMeshDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
