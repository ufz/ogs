/**
 * \file   CreateStructuredGridDialog.h
 * \author Karsten Rink
 * \date   2016-02-04
 * \brief  Definition of the CreateStructuredGridDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_CreateStructuredGrid.h"
#include <QDialog>

#include "MeshLib/Mesh.h"

/**
 * \brief A dialog window for managing general Data Explorer settings
 */
class CreateStructuredGridDialog : public QDialog, private Ui_CreateStructuredGrid
{
    Q_OBJECT

public:
    CreateStructuredGridDialog(QDialog* parent = nullptr);

private slots:
    void on_lineButton_toggled()  const;
    void on_triButton_toggled()   const { enable2dWidgets(); }
    void on_quadButton_toggled()  const { enable2dWidgets(); }
    void on_prismButton_toggled() const { enable3dWidgets(); }
    void on_hexButton_toggled()   const { enable3dWidgets(); }
    void on_meshExtentButton_toggled();
    void on_elemExtentButton_toggled();

    /// Instructions if the OK-Button has been pressed.
    void accept();

    /// Instructions if the Cancel-Button has been pressed.
    void reject() { this->done(QDialog::Rejected); };

private:
    void enable2dWidgets() const;
    void enable3dWidgets() const;
    void setValidators();

    /// Checks if all necessary inputs have been specified.
    bool inputIsEmpty() const;

signals:
    void meshAdded(MeshLib::Mesh* mesh);
};
