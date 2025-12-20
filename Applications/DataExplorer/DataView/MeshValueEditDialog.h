// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ui_MeshValueEdit.h"

#include <QDialog>

namespace MeshLib {
    class Mesh;
}

/**
 * \brief A dialog window for changing the MaterialID for mesh elements
 */
class MeshValueEditDialog : public QDialog, private Ui_MeshValueEdit
{
    Q_OBJECT

public:
    /// Constructor for creating a new FEM condition.
    explicit MeshValueEditDialog(MeshLib::Mesh* mesh,
                                 QDialog* parent = nullptr);

    ~MeshValueEditDialog() override;

private:
    MeshLib::Mesh* _mesh;

private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override;

    void on_replaceButton_toggled(bool isSelected);

signals:
    void valueEditFinished(MeshLib::Mesh*);
};
