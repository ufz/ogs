/**
 * \file
 * \author Lars Bilke
 * \date   2009-11-03
 * \brief  Definition of the StationTabWidget class.
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
#include "ui_StationTabWidgetBase.h"

/**
 * \brief Widget containing StationTreeView-objects.
 */
class StationTabWidget : public QWidget, public Ui_StationTabWidgetBase
{
    Q_OBJECT

public:
    StationTabWidget(QWidget* parent = nullptr);

private:

private slots:
    void enableSaveButton(bool enable) { this->saveStnPushButton->setEnabled(enable); };
    void enableRemoveButton(bool enable) { this->removeStnPushButton->setEnabled(enable); };
};
