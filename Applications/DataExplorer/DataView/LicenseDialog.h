/**
 * \file
 * \author Karsten Rink
 * \date   2012-11-30
 * \brief  Definition of the LicenseDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
    LicenseDialog(QDialog* parent = nullptr);
    ~LicenseDialog() override = default;
    ;

private:
    void setText();

private slots:
    void on_okPushButton_pressed();

};
