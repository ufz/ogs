/**
 * \file
 * \author Lars Bilke
 * \date   2009-11-03
 * \brief  Implementation of the MshTabWidget class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "MshTabWidget.h"

MshTabWidget::MshTabWidget( QWidget* parent /*= 0*/ )
    : QWidget(parent)
{
    setupUi(this);

    connect(this->addMeshPushButton, SIGNAL(clicked()), this->treeView, SLOT(addMesh()));
    connect(this->saveMeshPushButton, SIGNAL(clicked()), this->treeView, SLOT(writeToFile()));
    connect(this->removeMeshPushButton, SIGNAL(clicked()), this->treeView, SLOT(removeMesh()));
    connect(this->treeView, SIGNAL(enableSaveButton(bool)), this, SLOT(enableSaveButton(bool)));
    connect(this->treeView, SIGNAL(enableRemoveButton(bool)), this, SLOT(enableRemoveButton(bool)));
}
