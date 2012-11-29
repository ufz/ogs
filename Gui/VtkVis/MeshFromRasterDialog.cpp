/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file MeshFromRasterDialog.cpp
 *
 * Created on 2011-11-17 by Karsten Rink
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

