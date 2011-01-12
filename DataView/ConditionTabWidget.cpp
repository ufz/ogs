/**
 * \file ConditionTabWidget.cpp
 * 2010/12/14 KR Initial implementation
 *
 * Implementation of StationTabWidget
 */

// ** INCLUDES **
#include "ConditionTabWidget.h"
#include "TreeItem.h"
#include "ConditionModel.h"

ConditionTabWidget::ConditionTabWidget( QWidget* parent /*= 0*/ )
: QWidget(parent)
{
	setupUi(this);
}

void ConditionTabWidget::removeCondition()
{
	emit requestConditionRemoval(this->treeView->selectionModel()->currentIndex());
}

void ConditionTabWidget::removeAllConditions()
{
	ConditionModel* model = static_cast<ConditionModel*>(this->treeView->model());

	for (size_t j=0; j<3; j++)
	{
		QModelIndex parentIndex = model->index(j, 0, QModelIndex());
		int nChildren = model->getItem(parentIndex)->childCount();
		for (int i=nChildren; i>=0; i--)
			emit requestConditionRemoval(model->index(i, 0, parentIndex));
	}
}


void ConditionTabWidget::contextMenuEvent( QContextMenuEvent* event )
{
	/*
	QMenu menu;
	QAction* editMeshAction    = menu.addAction("Edit mesh...");
	QAction* saveMeshAction    = menu.addAction("Save mesh...");
	connect(editMeshAction, SIGNAL(triggered()), this, SLOT(openMshEditDialog()));
	connect(saveMeshAction, SIGNAL(triggered()), this, SLOT(writeMeshToFile()));
	menu.exec(event->globalPos());
	*/
}
