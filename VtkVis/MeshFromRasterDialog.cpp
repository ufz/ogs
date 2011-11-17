/**
 * \file MeshFromRasterDialog.cpp
 * 2011/11/17 KR Initial implementation
 */

#include "MeshFromRasterDialog.h"
#include "OGSError.h"

MeshFromRasterDialog::MeshFromRasterDialog(QDialog* parent)
: QDialog(parent)
{
	setupUi(this);

	this->elevationButton->setChecked(true);
	this->triButton->setChecked(true);
	this->mshNameEdit->setText("NewMesh");
}


MeshFromRasterDialog::~MeshFromRasterDialog()
{
}


void MeshFromRasterDialog::accept()
{
	emit setMeshParameters(this->mshNameEdit->text(),
						   this->triButton->isChecked(),
						   this->elevationButton->isChecked());
	this->done(QDialog::Accepted);
}

void MeshFromRasterDialog::reject()
{
	this->done(QDialog::Rejected);
}

