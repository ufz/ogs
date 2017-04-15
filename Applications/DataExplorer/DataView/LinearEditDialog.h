/**
 * \file
 * \author Karsten Rink
 * \date   2012-04-17
 * \brief  Definition of the LinearEditDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_LinearEdit.h"
#include <QDialog>

#include "Polyline.h"

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
    ~LinearEditDialog(void);

private:
    void setupDialog(const std::vector<std::size_t> &dis_nodes, const std::vector<double> &dis_values);

    const GeoLib::Polyline _line;

private slots:
    void on_comboBox_currentIndexChanged(int index);

    /// Instructions if the OK-Button has been pressed.
    void accept();

    /// Instructions if the Cancel-Button has been pressed.
    void reject();

signals:
    void transmitDisValues(std::vector< std::pair<std::size_t,double> >);
};
