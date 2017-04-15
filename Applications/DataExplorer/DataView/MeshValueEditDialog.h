/**
 * \file
 * \author Karsten Rink
 * \date   2013-03-27
 * \brief  Definition of the MeshValueEditDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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

    ~MeshValueEditDialog(void);

private:
    MeshLib::Mesh* _mesh;

private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept();

    /// Instructions if the Cancel-Button has been pressed.
    void reject();

    void on_replaceButton_toggled(bool isSelected);

signals:
    void valueEditFinished(MeshLib::Mesh*);
};
