// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "AddLayerToMeshDialog.h"

#include "Base/OGSError.h"
#include "Base/StrictDoubleValidator.h"

AddLayerToMeshDialog::AddLayerToMeshDialog(QDialog* parent) : QDialog(parent)
{
    setupUi(this);

    auto* thickness_validator = new StrictDoubleValidator(0, 1000000, 7, this);
    this->thicknessEdit->setValidator(thickness_validator);
}

void AddLayerToMeshDialog::accept()
{
    if (this->nameEdit->text().isEmpty())
    {
        OGSError::box("Please enter a name for the new Mesh.");
    }
    else if (this->thicknessEdit->text().isEmpty() ||
             this->thicknessEdit->text().toDouble() <= 0)
    {
        OGSError::box("Thickness needs to be larger 0");
    }
    else
    {
        this->done(QDialog::Accepted);
    }
}

void AddLayerToMeshDialog::reject()
{
    this->done(QDialog::Rejected);
}
