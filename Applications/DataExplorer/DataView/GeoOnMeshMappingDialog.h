// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    explicit GeoOnMeshMappingDialog(
        std::vector<std::unique_ptr<MeshLib::Mesh>> const& mesh_vec,
        QDialog* parent = nullptr);
    ~GeoOnMeshMappingDialog() override;

    std::string const& getNewGeoName() const { return _new_geo_name; };
    int getDataSetChoice() const;

private:
    std::string _new_geo_name;

private slots:
    void on_advancedMappingButton_toggled(bool isSelected) { this->geoNameEdit->setEnabled(isSelected); };

    void on_meshNameComboBox_currentIndexChanged(int idx);

    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override { this->done(QDialog::Rejected); };
};
