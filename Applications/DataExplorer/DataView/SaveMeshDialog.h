/**
 * \file
 * \author Karsten Rink
 * \date   2014-10-27
 * \brief  Definition of the SaveMeshDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_SaveMesh.h"
#include <QDialog>

namespace MeshLib {
    class Mesh;
}

/**
 * \brief A dialog window for managing properties for writing meshes to files.
 */
class SaveMeshDialog : public QDialog, private Ui_SaveMesh
{
    Q_OBJECT

public:
    explicit SaveMeshDialog(MeshLib::Mesh const& mesh,
                            QDialog* parent = nullptr);
    ~SaveMeshDialog() override = default;

private slots:
    /// Selection of path to save file
    void on_selectDirButton_clicked();

    void on_dataModeBox_currentIndexChanged(int index);

    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override { this->done(QDialog::Rejected); };

private:
    MeshLib::Mesh const& _mesh;

};
