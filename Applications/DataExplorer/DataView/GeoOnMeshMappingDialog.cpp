/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-07
 * \brief  Implementation of the GeoOnMeshMappingDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "GeoOnMeshMappingDialog.h"
#include "Mesh.h"


#include "OGSError.h"

GeoOnMeshMappingDialog::GeoOnMeshMappingDialog(
    std::vector<std::unique_ptr<MeshLib::Mesh>> const& mesh_vec,
    QDialog* parent)
: QDialog(parent), new_geo_name_("")
{
    setupUi(this);

    for (const auto& mesh : mesh_vec)
    {
        this->meshNameComboBox->addItem(
            QString::fromStdString(mesh->getName()));
    }
}

GeoOnMeshMappingDialog::~GeoOnMeshMappingDialog() = default;

int GeoOnMeshMappingDialog::getDataSetChoice() const
{
    return this->meshNameComboBox->currentIndex();
}

void GeoOnMeshMappingDialog::on_meshNameComboBox_currentIndexChanged(int idx)
{
    if (idx == 1)
    {
        this->normalMappingButton->setChecked(true);
    }

    bool is_enabled(idx != 1);
    this->normalMappingButton->setEnabled(is_enabled);
    this->advancedMappingButton->setEnabled(is_enabled);
    this->geoNameEdit->setEnabled(is_enabled && this->advancedMappingButton->isChecked());
}

void GeoOnMeshMappingDialog::accept()
{
    if (this->advancedMappingButton->isChecked())
    {
        new_geo_name_ = this->geoNameEdit->text().toStdString();
        if (new_geo_name_.empty())
        {
            OGSError::box("Please enter name for new geometry.");
            return;
        }
    }
    this->done(QDialog::Accepted);
}


