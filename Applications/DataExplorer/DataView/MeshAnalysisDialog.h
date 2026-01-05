// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#include "ui_MeshAnalysis.h"
#include <QDialog>

#include "MeshLib/Elements/ElementErrorCode.h"

namespace MeshLib {
    class Mesh;
}

/**
 * \brief A dialog window for calling mesh analysis methods
 */
class MeshAnalysisDialog : public QDialog, private Ui_MeshAnalysis
{
    Q_OBJECT

public:
    explicit MeshAnalysisDialog(
        std::vector<std::unique_ptr<MeshLib::Mesh>> const& mesh_vec,
        QDialog* parent = nullptr);
    ~MeshAnalysisDialog() override;

private:
    /// Prepares the output for the node message window
    void nodesMsgOutput(std::vector<std::size_t> const& node_ids, std::vector<std::size_t> const& collapsibleNodeIds);

    /// Prepares the output for the node message window
    void elementsMsgOutput(const std::vector<ElementErrorCode> &error_codes);

    std::vector<std::unique_ptr<MeshLib::Mesh>> const& _mesh_vec;

private slots:
    /// Starts the analysis
    void on_startButton_pressed();

    /// Closes the dialog
    void on_closeButton_pressed() { this->close(); }
};
