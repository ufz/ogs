// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    explicit DataExplorerSettingsDialog(QDialog* parent = nullptr);
    ~DataExplorerSettingsDialog() override;

private slots:
    void on_gmshPathButton_clicked();

    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override { this->done(QDialog::Rejected); };
};
