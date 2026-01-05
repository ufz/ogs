// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ui_AddLayerToMesh.h"

#include <QDialog>
#include <QLineEdit>

/**
 * \brief A dialog window for adding a layer to the top or bottom of a mesh
 */
class AddLayerToMeshDialog : public QDialog, private Ui_AddLayerToMesh
{
    Q_OBJECT

public:
    explicit AddLayerToMeshDialog(QDialog* parent = nullptr);

    /// Returns if the top layer button is selected (if false, bottom is selected).
    bool isTopLayer() const { return this->topButton->isChecked(); };

    /// Returns the thickness of the new layer.
    double getThickness() const { return this->thicknessEdit->text().toDouble(); };

    /// Returns the name of the new mesh.
    std::string getName() const { return this->nameEdit->text().toStdString(); };

private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override;
};
