/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
