// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

// ** INCLUDES **
#include "ui_GeoTabWidgetBase.h"

/**
 * \brief Widget containing GeoTreeView-objects.
 */
class GeoTabWidget : public QWidget, public Ui_GeoTabWidgetBase
{
    Q_OBJECT

public:
    explicit GeoTabWidget(QWidget* parent = nullptr);

private:

private slots:
    void enableSaveButton(bool enable) { this->saveGeoPushButton->setEnabled(enable); };
    void enableRemoveButton(bool enable) { this->removeGeoPushButton->setEnabled(enable); };

};
