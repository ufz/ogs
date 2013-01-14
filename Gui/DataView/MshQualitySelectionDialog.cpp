/**
 * \file
 * \author Karsten Rink
 * \date   2011-03-16
 * \brief  Implementation of the MshQualitySelectionDialog class.
 *
 * \copyright
 * Copyright (c)  2013, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MshQualitySelectionDialog.h"
#include "VtkMeshSource.h"

/// Constructor
MshQualitySelectionDialog::MshQualitySelectionDialog(VtkMeshSource* msh, QDialog* parent)
	: QDialog(parent), _msh(msh)
{
	setupUi(this);
	this->choiceEdges->toggle();
}

MshQualitySelectionDialog::~MshQualitySelectionDialog()
{
}

/// Instructions if the OK-Button has been pressed.
void MshQualitySelectionDialog::accept()
{
	MshQualityType::type t;
	if (this->choiceEdges->isChecked())
		t = MshQualityType::EDGERATIO;
	else if (this->choiceArea->isChecked())
		t = MshQualityType::AREA;
	else if (this->choiceVolume->isChecked())
		t = MshQualityType::VOLUME;
	else if (this->choiceAngles->isChecked())
		t = MshQualityType::EQUIANGLESKEW;
	else
		t = MshQualityType::INVALID;

	emit measureSelected(_msh, t);
	this->done(QDialog::Accepted);
}

/// Instructions if the Cancel-Button has been pressed.
void MshQualitySelectionDialog::reject()
{
	this->done(QDialog::Rejected);
}
