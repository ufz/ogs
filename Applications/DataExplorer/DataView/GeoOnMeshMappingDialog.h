/**
 * \file
 * \author Karsten Rink
 * \date   2013-10-07
 * \brief  Definition of the GeoOnMeshMappingDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_GeoOnMeshMapping.h"
#include <QDialog>

#include "MeshGeoToolsLib/GeoMapper.h"

namespace MeshLib {
    class Mesh;
}

/**
 * \brief A dialog window for creating DIRECT boundary conditions from raster files
 */
class GeoOnMeshMappingDialog : public QDialog, private Ui_GeoOnMeshMapping
{
    Q_OBJECT

public:
    GeoOnMeshMappingDialog(
        std::vector<std::unique_ptr<MeshLib::Mesh>> const& mesh_vec,
        QDialog* parent = 0);
    ~GeoOnMeshMappingDialog(void);

    std::string const& getNewGeoName() const { return _new_geo_name; };
    int getDataSetChoice() const;

private:
    std::string _new_geo_name;

private slots:
    void on_advancedMappingButton_toggled(bool isSelected) { this->geoNameEdit->setEnabled(isSelected); };

    void on_meshNameComboBox_currentIndexChanged(int idx);

    /// Instructions if the OK-Button has been pressed.
    void accept();

    /// Instructions if the Cancel-Button has been pressed.
    void reject() { this->done(QDialog::Rejected); };

};
