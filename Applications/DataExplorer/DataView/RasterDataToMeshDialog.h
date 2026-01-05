// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ui_RasterDataToMesh.h"

#include <QDialog>
#include <QLineEdit>
#include <QString>

/**
 * \brief A dialog window for transferring raster data onto a mesh.
 */
class RasterDataToMeshDialog : public QDialog, private Ui_RasterDataToMesh
{
    Q_OBJECT

public:
    explicit RasterDataToMeshDialog(std::string const& mesh_name,
                                    QDialog* parent = nullptr);

    bool createNodeArray() const { return this->nodeButton->isChecked(); }
    bool createElemArray() const { return this->elemButton->isChecked(); }
    std::string getArrayName() const { return this->arrayNameEdit->text().toStdString(); }
    std::string getMeshName() const { return this->meshNameEdit->text().toStdString(); }
    std::string getRasterPath() const { return this->rasterPathEdit->text().toStdString(); }
    double getNoDataReplacement() const { return this->noDataValueEdit->text().toDouble(); }

private slots:
    void on_rasterSelectButton_pressed();

    void accept() override;
    void reject() override { this->done(QDialog::Rejected); }
};
