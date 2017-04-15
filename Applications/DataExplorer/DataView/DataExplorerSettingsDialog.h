/**
 * \file   DataExplorerSettingsDialog.h
 * \author Karsten Rink
 * \date   2014-02-05
 * \brief  Definition of the DataExplorerSettingsDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_DataExplorerSettings.h"
#include <QDialog>

/**
 * \brief A dialog window for managing general Data Explorer settings
 */
class DataExplorerSettingsDialog : public QDialog, private Ui_DataExplorerSettings
{
    Q_OBJECT

public:
    DataExplorerSettingsDialog(QDialog* parent = nullptr);
    ~DataExplorerSettingsDialog(void) override;

private slots:
    void on_gmshPathButton_clicked();

    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override { this->done(QDialog::Rejected); };
};
