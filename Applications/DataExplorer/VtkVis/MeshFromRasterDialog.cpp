/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-17
 * \brief  Implementation of the MeshFromRasterDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshFromRasterDialog.h"
#include "OGSError.h"
#include "MeshGenerators/VtkMeshConverter.h"

MeshFromRasterDialog::MeshFromRasterDialog(QDialog* parent)
: QDialog(parent), _mesh_name("mesh"), _array_name("MaterialIDs")
{
    setupUi(this);

    this->elevationButton->setChecked(true);
    this->triButton->setChecked(true);
    this->mshNameEdit->setText("RasterDataMesh");
}

MeshFromRasterDialog::~MeshFromRasterDialog() = default;

void MeshFromRasterDialog::on_elevationButton_toggled(bool isChecked)
{
    if (isChecked)
    {
        if (this->prismButton->isChecked())
            this->triButton->setChecked(true);
        if (this->hexButton->isChecked())
            this->quadButton->setChecked(true);
    }

    this->prismButton->setEnabled(!isChecked);
    this->hexButton->setEnabled(!isChecked);
}

void MeshFromRasterDialog::accept()
{
    if (this->mshNameEdit->text().isEmpty())
    {
        OGSError::box("Please specify a name for the resulting mesh.");
        return;
    }
    _mesh_name = this->mshNameEdit->text().toStdString();

    _intensity_selection = MeshLib::UseIntensityAs::ELEVATION;
    if (this->materialButton->isChecked()) _intensity_selection  = MeshLib::UseIntensityAs::MATERIALS;
    else if (this->otherButton->isChecked()) _intensity_selection  = MeshLib::UseIntensityAs::DATAVECTOR;
    else if (this->ignoreButton->isChecked()) _intensity_selection  = MeshLib::UseIntensityAs::NONE;

    if (_intensity_selection == MeshLib::UseIntensityAs::DATAVECTOR)
    {
        if (this->arrayNameEdit->text().isEmpty())
        {
            OGSError::box("Please specify a name for the data vector.");
            return;
        }

        _array_name = this->arrayNameEdit->text().toStdString();
    }
    _element_selection = MeshLib::MeshElemType::TRIANGLE;
    if (this->quadButton->isChecked()) _element_selection = MeshLib::MeshElemType::QUAD;
    else if (this->prismButton->isChecked()) _element_selection = MeshLib::MeshElemType::PRISM;
    else if (this->hexButton->isChecked()) _element_selection = MeshLib::MeshElemType::HEXAHEDRON;

    this->done(QDialog::Accepted);
}

void MeshFromRasterDialog::reject()
{
    this->done(QDialog::Rejected);
}

