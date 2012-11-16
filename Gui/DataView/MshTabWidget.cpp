/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file MshTabWidget.cpp
 *
 * Created on 2009-11-03 by Lars Bilke
 */

// ** INCLUDES **
#include "MshTabWidget.h"

MshTabWidget::MshTabWidget( QWidget* parent /*= 0*/ )
	: QWidget(parent)
{
	setupUi(this);

	connect(this->addMeshPushButton, SIGNAL(clicked()), this->treeView, SLOT(addMeshAction()));
	connect(this->saveMeshPushButton, SIGNAL(clicked()), this->treeView, SLOT(writeMeshToFile()));
	connect(this->removeMeshPushButton, SIGNAL(clicked()), this->treeView, SLOT(removeMesh()));
	//connect(this->clearAllPushButton, SIGNAL(clicked()), this->treeView, SLOT(removeAllMeshes()));
}
