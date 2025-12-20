// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    explicit StationTabWidget(QWidget* parent = nullptr);

private:

private slots:
    void enableSaveButton(bool enable) { this->saveStnPushButton->setEnabled(enable); };
    void enableRemoveButton(bool enable) { this->removeStnPushButton->setEnabled(enable); };
};
