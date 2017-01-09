/**
 * \file
 * \author Karsten Rink
 * \date   2012-01-04
 * \brief  Definition of the CondFromRasterDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CONDFROMRASTERDIALOG_H
#define CONDFROMRASTERDIALOG_H

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
    CondFromRasterDialog(const std::vector<MeshLib::Mesh*> &msh_vec, QDialog* parent = 0);
    ~CondFromRasterDialog(void);

private:
    const std::vector<MeshLib::Mesh*> _msh_vec;
    StrictDoubleValidator* _scale_validator;

private slots:
    void on_integrateButton_toggled(bool isSelected);
    void on_selectButton_pressed();

    /// Instructions if the OK-Button has been pressed.
    void accept();

    /// Instructions if the Cancel-Button has been pressed.
    void reject();

signals:
    void directNodesWritten(std::string);
    void transmitDisValues(std::vector< std::pair<std::size_t,double> >);
};

#endif //CONDFROMRASTERDIALOG_H
