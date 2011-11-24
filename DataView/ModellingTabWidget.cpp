/**
 * \file ModellingTabWidget.cpp
 * 2010/12/14 KR Initial implementation
 *
 * Implementation of StationTabWidget
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
