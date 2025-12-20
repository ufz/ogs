// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ui_LinearEdit.h"
#include <QDialog>

#include "GeoLib/Polyline.h"

/**
 * \brief A dialog window for creating linear boundary conditions on polylines
 */
class LinearEditDialog : public QDialog, private Ui_LinearEdit
{
    Q_OBJECT

public:
    LinearEditDialog(const GeoLib::Polyline& line,
                     const std::vector<std::size_t>& dis_nodes,
                     const std::vector<double>& dis_values,
                     QDialog* parent = nullptr);
    ~LinearEditDialog() override;

private:
    void setupDialog(const std::vector<std::size_t> &dis_nodes, const std::vector<double> &dis_values);

    const GeoLib::Polyline _line;

private slots:
    void on_comboBox_currentIndexChanged(int index);

    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override;

signals:
    void transmitDisValues(std::vector< std::pair<std::size_t,double> >);
};
