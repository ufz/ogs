/**
 * \file
 * \author Karsten Rink
 * \date   2011-02-07
 * \brief  Implementation of the GeoTabWidget class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "GeoTabWidget.h"

GeoTabWidget::GeoTabWidget( QWidget* parent /*= 0*/ )
    : QWidget(parent)
{
    setupUi(this);

    connect(this->openGeoPushButton, SIGNAL(clicked()), this->treeView, SLOT(addGeometry()));
    connect(this->saveGeoPushButton, SIGNAL(clicked()), this->treeView, SLOT(writeToFile()));
    connect(this->removeGeoPushButton, SIGNAL(clicked()), this->treeView, SLOT(removeGeometry()));
    connect(this->treeView, SIGNAL(enableSaveButton(bool)), this, SLOT(enableSaveButton(bool)));
    connect(this->treeView, SIGNAL(enableRemoveButton(bool)), this, SLOT(enableRemoveButton(bool)));
}
