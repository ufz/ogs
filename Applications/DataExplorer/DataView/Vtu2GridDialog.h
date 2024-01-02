/**
 * \file
 * \date   2023-04-26
 * \brief  Definition of the Vtu2GridDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QDialog>
#include <QStringListModel>
#include <memory>

#include "MeshLib/Elements/ElementErrorCode.h"
#include "ui_Vtu2Grid.h"

class MeshModel;

/*
 * \brief A dialog window for calling methods to create a 3D Voxelgrid from
 * multiple 2D vtu meshes
 */
class Vtu2GridDialog : public QDialog, private Ui_Vtu2Grid
{
    Q_OBJECT

public:
    explicit Vtu2GridDialog(MeshModel& mesh_model, QDialog* parent = nullptr);

private:
    MeshModel& _mesh_model;
    QStringListModel _allMeshes;

private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept() override;
    /// Instructions if the Cancel-Button has been pressed.
    void reject() override { this->done(QDialog::Rejected); };
    /// As the x/y/z input changes an estimation of the expected Voxel is given.
    void updateExpectedVoxel();
    void on_xlineEdit_textChanged();
    void on_ylineEdit_textChanged();
    void on_zlineEdit_textChanged();
};