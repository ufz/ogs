/**
 * \file
 * \author Karsten Rink
 * \date   2011-03-16
 * \brief  Implementation of the MshQualitySelectionDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshQualitySelectionDialog.h"
#include "VtkMeshSource.h"

/// Constructor
MshQualitySelectionDialog::MshQualitySelectionDialog(QDialog* parent)
: QDialog(parent), _metric (MeshQualityType::EDGERATIO)
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
	if (this->choiceEdges->isChecked())
		_metric = MeshQualityType::EDGERATIO;
	else if (this->choiceArea->isChecked())
		_metric = MeshQualityType::AREA;
	else if (this->choiceVolume->isChecked())
		_metric = MeshQualityType::VOLUME;
	else if (this->choiceAngles->isChecked())
		_metric = MeshQualityType::EQUIANGLESKEW;
	else if (this->choiceRadius->isChecked())
		_metric = MeshQualityType::RADIUSEDGERATIO;
	else
		_metric = MeshQualityType::INVALID;

	this->done(QDialog::Accepted);
}

/// Instructions if the Cancel-Button has been pressed.
void MshQualitySelectionDialog::reject()
{
	this->done(QDialog::Rejected);
}
