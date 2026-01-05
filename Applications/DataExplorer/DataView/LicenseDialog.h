// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ui_License.h"
#include <QDialog>

/**
 * \brief A dialog window displaying the OGS license information
 */
class LicenseDialog : public QDialog, private Ui_License
{
    Q_OBJECT

public:
    explicit LicenseDialog(QDialog* parent = nullptr);
    ~LicenseDialog() override = default;
    ;

private:
    void setText();

private slots:
    void on_okPushButton_pressed();

};
