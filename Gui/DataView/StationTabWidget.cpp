/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file StationTabWidget.cpp
 *
 * Created on 2009-11-03 by Lars Bilke
 */

// ** INCLUDES **
#include "StationTabWidget.h"

StationTabWidget::StationTabWidget( QWidget* parent /*= 0*/ )
	: QWidget(parent)
{
	setupUi(this);

	connect(this->openStnPushButton, SIGNAL(clicked()), this->treeView, SLOT(addStationList()));
	connect(this->saveStnPushButton, SIGNAL(clicked()), this->treeView, SLOT(writeToFile()));
	connect(this->removeStnPushButton, SIGNAL(clicked()), this->treeView, SLOT(removeStationList()));
	connect(this->treeView, SIGNAL(enableSaveButton(bool)), this, SLOT(enableSaveButton(bool)));
	connect(this->treeView, SIGNAL(enableRemoveButton(bool)), this, SLOT(enableRemoveButton(bool)));
}
