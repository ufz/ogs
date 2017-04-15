/**
 * \file   AddLayerToMeshDialog.cpp
 * \author Karsten Rink
 * \date   2016-01-18
 * \brief  Implementation of the AddLayerToMeshDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AddLayerToMeshDialog.h"
#include "OGSError.h"
#include "StrictDoubleValidator.h"


AddLayerToMeshDialog::AddLayerToMeshDialog(QDialog* parent)
: QDialog(parent)
{
    setupUi(this);

    auto* thickness_validator = new StrictDoubleValidator(0, 1000000, 7, this);
    this->thicknessEdit->setValidator (thickness_validator);
}

void AddLayerToMeshDialog::accept()
{
    if (this->nameEdit->text().isEmpty())
        OGSError::box("Please enter a name for the new Mesh.");
    else if (this->thicknessEdit->text().isEmpty() ||
        this->thicknessEdit->text().toDouble() <= 0)
        OGSError::box("Thickness needs to be larger 0");
    else
        this->done(QDialog::Accepted);
}

void AddLayerToMeshDialog::reject()
{
    this->done(QDialog::Rejected);
}

