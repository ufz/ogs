/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-17
 * \brief  Implementation of the MeshFromRasterDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshFromRasterDialog.h"
#include "OGSError.h"
#include "MeshGenerators/VtkMeshConverter.h"

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
	MeshLib::UseIntensityAs i_type(MeshLib::UseIntensityAs::ELEVATION);
	if (this->materialButton->isChecked()) i_type = MeshLib::UseIntensityAs::DATAVECTOR;
	else if (this->ignoreButton->isChecked()) i_type = MeshLib::UseIntensityAs::NONE;

	MeshLib::MeshElemType e_type(MeshLib::MeshElemType::TRIANGLE);
	if (this->quadButton->isChecked()) e_type = MeshLib::MeshElemType::QUAD;
	else if (this->hexButton->isChecked()) e_type = MeshLib::MeshElemType::HEXAHEDRON;

	emit setMeshParameters(this->mshNameEdit->text(), e_type, i_type);
	this->done(QDialog::Accepted);
}

void MeshFromRasterDialog::reject()
{
	this->done(QDialog::Rejected);
}

