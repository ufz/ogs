/**
 * \file
 * \date   2023-04-14
 * \brief  Definition of the TranslateDataDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <QDialog>
#include <memory>

#include "MeshLib/Elements/ElementErrorCode.h"
#include "MeshLib/MeshEditing/moveMeshNodes.h"
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
    ~TranslateDataDialog() override;

private:
    MeshModel* _mesh_model;
    GEOModels* _geo_models;

private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept() override;
    /// Instructions if the Cancel-Button has been pressed.
    void reject() override { this->done(QDialog::Rejected); };

    /*
    private:
        /// Prepares the output for the node message window
        void nodesMsgOutput(std::vector<std::size_t> const& node_ids,
    std::vector<std::size_t> const& collapsibleNodeIds);

        /// Prepares the output for the node message window
        void elementsMsgOutput(const std::vector<ElementErrorCode>
    &error_codes);

        std::vector<std::unique_ptr<MeshLib::Mesh>> const& _mesh_vec;

    private slots:
        /// Starts the analysis
        void on_startButton_pressed();

        /// Closes the dialog
        void on_closeButton_pressed() { this->close(); }
    */
signals:
    void meshAdded(MeshLib::Mesh const* mesh);
};
