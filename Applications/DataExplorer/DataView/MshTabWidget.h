/**
 * \file
 * \author Lars Bilke
 * \date   2009-11-03
 * \brief  Definition of the MshTabWidget class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

// ** INCLUDES **
#include "ui_MshTabWidgetBase.h"

/**
 * \brief Widget for data views of meshes.
 */
class MshTabWidget : public QWidget, public Ui_MshTabWidgetBase
{
    Q_OBJECT

public:
    MshTabWidget(QWidget* parent = nullptr);

private slots:
    void enableSaveButton(bool enable) { this->saveMeshPushButton->setEnabled(enable); };
    void enableRemoveButton(bool enable) { this->removeMeshPushButton->setEnabled(enable); };
};
