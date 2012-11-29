/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file ModellingTabWidget.cpp
 *
 * Created on 2010-12-14 by Karsten Rink
 */

// ** INCLUDES **
#include "ProcessModel.h"
#include "ModellingTabWidget.h"

ModellingTabWidget::ModellingTabWidget( QWidget* parent /*= 0*/ )
	: QWidget(parent)
{
	setupUi(this);
}

void ModellingTabWidget::on_addProcessButton_pressed()
{
	emit requestNewProcess();
}

void ModellingTabWidget::on_deleteAllButton_pressed()
{
	static_cast<ProcessModel*>(this->treeView->model())->removeAllProcesses();
}
