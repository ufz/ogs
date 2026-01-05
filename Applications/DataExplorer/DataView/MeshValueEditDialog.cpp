// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MeshValueEditDialog.h"

#include "Base/OGSError.h"
#include "MeshToolsLib/MeshEditing/ElementValueModification.h"

MeshValueEditDialog::MeshValueEditDialog(MeshLib::Mesh* mesh, QDialog* parent)
    : QDialog(parent), _mesh(mesh)
{
    setupUi(this);
    this->edit_old_value->setEnabled(false);
    this->edit_new_value->setEnabled(false);
    this->replaceCheckBox->setEnabled(false);
}

MeshValueEditDialog::~MeshValueEditDialog(void) = default;

void MeshValueEditDialog::accept()
{
    if (this->condenseButton->isChecked())
    {
        MeshToolsLib::ElementValueModification::condense(*_mesh);
    }
    else
    {
        if (this->edit_old_value->text().isEmpty())
        {
            OGSError::box("Please input which material you want to replace.");
            return;
        }
        unsigned old_value =
            static_cast<unsigned>(this->edit_old_value->text().toInt());
        if (this->edit_new_value->text().isEmpty())
        {
            OGSError::box("Please input the new material to replace group " +
                          this->edit_old_value->text() + ".");
            return;
        }
        unsigned new_value =
            static_cast<unsigned>(this->edit_new_value->text().toInt());
        bool do_not_replace = this->replaceCheckBox->isChecked();
        bool result = MeshToolsLib::ElementValueModification::replace(
            *_mesh, old_value, new_value, !do_not_replace);
        if (!result && do_not_replace)
        {
            OGSError::box("The new material group already exists.");
            return;
        }
    }

    emit valueEditFinished(_mesh);
    this->done(QDialog::Accepted);
}

void MeshValueEditDialog::reject()
{
    this->done(QDialog::Rejected);
}

void MeshValueEditDialog::on_replaceButton_toggled(bool isSelected)
{
    this->edit_old_value->setEnabled(isSelected);
    this->edit_new_value->setEnabled(isSelected);
    this->replaceCheckBox->setEnabled(isSelected);
}
