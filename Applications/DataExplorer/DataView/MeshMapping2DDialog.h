/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_MeshMapping2D.h"

#include <QDialog>
#include <QLineEdit>
#include <QString>

/**
 * \brief A dialog window for mapping a 2d mesh based on a raster file.
 */
class MeshMapping2DDialog : public QDialog, private Ui_MeshMapping2D
{
    Q_OBJECT

public:
    explicit MeshMapping2DDialog(QDialog* parent = nullptr);

    bool useRasterMapping() const { return this->rasterValueButton->isChecked(); }
    bool useStaticMapping() const { return this->staticValueButton->isChecked(); }
    std::string getRasterPath() const { return this->rasterPathEdit->text().toStdString(); }
    double getNoDataReplacement() const { return this->noDataValueEdit->text().toDouble(); }
    double getStaticValue() const { return this->staticValueEdit->text().toDouble(); }
    std::string getNewMeshName() const { return this->newNameEdit->text().toStdString(); }

private slots:
    void on_rasterValueButton_toggled(bool isChecked);
    void on_rasterSelectButton_pressed();

    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override { this->done(QDialog::Rejected); }
};
