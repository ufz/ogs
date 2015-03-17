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

#include "MshQualitySelectionDialog.h"
#include "vtkUnstructuredGridAlgorithm.h"

/// Constructor
MshQualitySelectionDialog::MshQualitySelectionDialog(InSituLib::VtkMappedMeshSource* mesh, QDialog* parent)
	: QDialog(parent), _mesh(mesh)
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
    MeshQualityType t;
    if (this->choiceEdges->isChecked())
        t = MeshQualityType::EDGERATIO;
    else if (this->choiceArea->isChecked())
        t = MeshQualityType::AREA;
    else if (this->choiceVolume->isChecked())
        t = MeshQualityType::VOLUME;
    else if (this->choiceAngles->isChecked())
        t = MeshQualityType::EQUIANGLESKEW;
    else if (this->choiceRadius->isChecked())
        t = MeshQualityType::RADIUSEDGERATIO;
    else
        t = MeshQualityType::INVALID;

    emit measureSelected(_mesh, t);
    this->done(QDialog::Accepted);
}

/// Instructions if the Cancel-Button has been pressed.
void MshQualitySelectionDialog::reject()
{
	this->done(QDialog::Rejected);
}
