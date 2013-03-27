/**
 * \file
 * \author Karsten Rink
 * \date   2013-03-27
 * \brief  Implementation of the MeshValueEditDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshValueEditDialog.h"
#include "MshEditor.h"
#include "OGSError.h"


MeshValueEditDialog::MeshValueEditDialog(MeshLib::Mesh* mesh, QDialog* parent)
	: QDialog(parent), _mesh(mesh)
{
	setupUi(this);
	this->edit_old_value->setEnabled(false);
	this->edit_new_value->setEnabled(false);
	this->replaceCheckBox->setEnabled(false);
}

MeshValueEditDialog::~MeshValueEditDialog(void)
{
}

void MeshValueEditDialog::accept()
{
	if (this->condenseButton->isChecked())
		MeshLib::MshEditor::condenseElementValues(*_mesh);
	else
	{
		if (this->edit_old_value->text().isEmpty())
		{
			OGSError::box("Please input which material you want to replace.");
			return;
		}
		unsigned old_value = static_cast<unsigned>(atoi(this->edit_old_value->text().toStdString().c_str()));
		if (this->edit_new_value->text().isEmpty())
		{
			OGSError::box("Please input the new material to replace group " + this->edit_old_value->text() + ".");
			return;
		}
		unsigned new_value = static_cast<unsigned>(atoi(this->edit_new_value->text().toStdString().c_str()));
		bool do_not_replace = this->replaceCheckBox->isChecked();
		bool result = MeshLib::MshEditor::replaceElementValue(*_mesh, old_value, new_value, !do_not_replace);
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


