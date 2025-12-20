// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

// ** INCLUDES **
#include "ui_MeshTabWidgetBase.h"

/**
 * \brief Widget for data views of meshes.
 */
class MeshTabWidget : public QWidget, public Ui_MeshTabWidgetBase
{
    Q_OBJECT

public:
    explicit MeshTabWidget(QWidget* parent = nullptr);

private slots:
    void enableSaveButton(bool enable) { this->saveMeshPushButton->setEnabled(enable); };
    void enableRemoveButton(bool enable) { this->removeMeshPushButton->setEnabled(enable); };
};
