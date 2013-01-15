/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-17
 * \brief  Implementation of the MeshFromRasterDialog class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshFromRasterDialog.h"
#include "OGSError.h"

MeshFromRasterDialog::MeshFromRasterDialog(QDialog* parent)
: QDialog(parent)
{
	setupUi(this);

	this->elevationButton->setChecked(true);
	this->triButton->setChecked(true);
	this->mshNameEdit->setText("RasterDataMesh");
}


MeshFromRasterDialog::~MeshFromRasterDialog()
{
}


void MeshFromRasterDialog::accept()
{
	UseIntensityAs::type i_type(UseIntensityAs::ELEVATION);
	if (this->materialButton->isChecked()) i_type = UseIntensityAs::MATERIAL;
	else if (this->ignoreButton->isChecked()) i_type = UseIntensityAs::NONE;

	MshElemType::type e_type(MshElemType::TRIANGLE);
	if (this->quadButton->isChecked()) e_type = MshElemType::QUAD;
	else if (this->hexButton->isChecked()) e_type = MshElemType::HEXAHEDRON;

	emit setMeshParameters(this->mshNameEdit->text(), e_type, i_type);
	this->done(QDialog::Accepted);
}

void MeshFromRasterDialog::reject()
{
	this->done(QDialog::Rejected);
}

