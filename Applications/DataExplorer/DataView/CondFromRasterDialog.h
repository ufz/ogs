// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ui_CondFromRaster.h"
#include <QDialog>

namespace MeshLib {
    class Mesh;
}

class StrictDoubleValidator;

/**
 * \brief A dialog window for creating DIRECT boundary conditions from raster files
 */
class CondFromRasterDialog : public QDialog, private Ui_CondFromRaster
{
    Q_OBJECT

public:
    explicit CondFromRasterDialog(std::vector<MeshLib::Mesh*> msh_vec,
                                  QDialog* parent = nullptr);
    ~CondFromRasterDialog() override;

private:
    const std::vector<MeshLib::Mesh*> _msh_vec;
    StrictDoubleValidator* _scale_validator;

private slots:
    void on_integrateButton_toggled(bool isSelected);
    void on_selectButton_pressed();

    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override;

signals:
    void directNodesWritten(std::string);
    void transmitDisValues(std::vector< std::pair<std::size_t,double> >);
};
