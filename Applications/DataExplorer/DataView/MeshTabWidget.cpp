// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// ** INCLUDES **
#include "MeshTabWidget.h"

MeshTabWidget::MeshTabWidget(QWidget* parent /*= 0*/) : QWidget(parent)
{
    setupUi(this);

    connect(this->addMeshPushButton, SIGNAL(clicked()), this->treeView,
            SLOT(addMesh()));
    connect(this->saveMeshPushButton, SIGNAL(clicked()), this->treeView,
            SLOT(writeToFile()));
    connect(this->removeMeshPushButton, SIGNAL(clicked()), this->treeView,
            SLOT(removeMesh()));
    connect(this->treeView, SIGNAL(enableSaveButton(bool)), this,
            SLOT(enableSaveButton(bool)));
    connect(this->treeView, SIGNAL(enableRemoveButton(bool)), this,
            SLOT(enableRemoveButton(bool)));
}
