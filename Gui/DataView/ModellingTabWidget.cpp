/**
 * \file
 * \author Karsten Rink
 * \date   2010-12-14
 * \brief  Implementation of the ModellingTabWidget class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
