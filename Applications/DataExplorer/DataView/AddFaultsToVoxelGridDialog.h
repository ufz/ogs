/**
 * \file
 * \date   2023-04-26
 * \brief  Definition of the AddFaultsToVoxelGridDialog class.
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
#include "ui_AddFaultsToVoxelGrid.h"

class MeshModel;

/*
 * \brief A dialog window for calling methods to include one or multiple
 * faults in a 3D voxelgrid.
 */
class AddFaultsToVoxelGridDialog : public QDialog,
                                   private Ui_AddFaultsToVoxelGrid
{
    Q_OBJECT

public:
    explicit AddFaultsToVoxelGridDialog(MeshModel& mesh_model,
                                        QDialog* parent = nullptr);

private:
    MeshModel& _mesh_model;
    QStringListModel _voxelGrids;
    QStringListModel _meshes2D;
    QStringListModel _selFaults;

private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept() override;
    /// Instructions if the Cancel-Button has been pressed.
    void reject() override { this->done(QDialog::Rejected); };
    /// Instructions if the ">>-button" has been pressed.
    void on_selectDataButton_pressed();
    /// Instructions if the "<<-button" has been pressed.
    void on_deselectDataButton_pressed();
};