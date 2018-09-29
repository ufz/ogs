/**
 * \file
 * \author Karsten Rink
 * \date   2013-03-27
 * \brief  Definition of the MeshValueEditDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
    MeshValueEditDialog(MeshLib::Mesh* mesh, QDialog* parent = nullptr);

    ~MeshValueEditDialog(void) override;

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
