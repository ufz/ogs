// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// ** INCLUDES **
#include "StationTabWidget.h"

StationTabWidget::StationTabWidget(QWidget* parent /*= 0*/) : QWidget(parent)
{
    setupUi(this);

    connect(this->openStnPushButton, SIGNAL(clicked()), this->treeView,
            SLOT(addStationList()));
    connect(this->saveStnPushButton, SIGNAL(clicked()), this->treeView,
            SLOT(writeToFile()));
    connect(this->removeStnPushButton, SIGNAL(clicked()), this->treeView,
            SLOT(removeStationList()));
    connect(this->treeView, SIGNAL(enableSaveButton(bool)), this,
            SLOT(enableSaveButton(bool)));
    connect(this->treeView, SIGNAL(enableRemoveButton(bool)), this,
            SLOT(enableRemoveButton(bool)));
}
