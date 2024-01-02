/**
 * \file
 * \date   2023-04-14
 * \brief  Definition of the TranslateDataDialog class.
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
#include "MeshToolsLib/MeshEditing/moveMeshNodes.h"
#include "ui_TranslateData.h"

namespace MeshLib
{
class Mesh;
}

class MeshModel;

class GEOModels;

/**
 * \brief A dialog window for calling translation methods
 */
class TranslateDataDialog : public QDialog, private Ui_TranslateData
{
    Q_OBJECT

public:
    explicit TranslateDataDialog(MeshModel* mesh_model,
                                 GEOModels* geo_models,
                                 QDialog* parent = nullptr);

private:
    void moveGeometry(Eigen::Vector3d const& displacement,
                      std::string const& name);
    void moveMesh(Eigen::Vector3d const& displacement, std::string const& name);
    MeshModel* _mesh_model;
    GEOModels* _geo_models;
    QStringListModel _allData;
    QStringListModel _selData;

private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept() override;
    /// Instructions if the Cancel-Button has been pressed.
    void reject() override { this->done(QDialog::Rejected); };
    /// Instructions if the ">>-button" has been pressed.
    void on_selectDataButton_pressed();
    /// Instructions if the "<<-button" has been pressed.
    void on_deselectDataButton_pressed();

signals:
    void meshAdded(MeshLib::Mesh const* mesh);
};
