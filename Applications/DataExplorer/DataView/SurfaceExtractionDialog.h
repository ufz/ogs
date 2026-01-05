// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <QDialog>

#include "ui_SurfaceExtraction.h"

namespace MeshLib {
    class Mesh;
}

/**
 * \brief A dialog window for managing properties for writing meshes to files.
 */
class SurfaceExtractionDialog : public QDialog, private Ui_SurfaceExtraction
{
    Q_OBJECT

public:
    explicit SurfaceExtractionDialog(QDialog* parent = nullptr);
    ~SurfaceExtractionDialog() override = default;

    int getTolerance() const { return _tolerance; }
    Eigen::Vector3d const& getNormal() const { return _dir; }

private slots:
    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override { this->done(QDialog::Rejected); };

private:
    int _tolerance{90};
    Eigen::Vector3d _dir{0, 0, -1};
};
