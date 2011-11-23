/**
 * \file ProcessView.cpp
 * 2010/12/13 KR Initial implementation
 */

#include <QMenu>

#include "ProcessItem.h"
#include "CondObjectListItem.h"
#include "ProcessModel.h"
#include "ProcessView.h"

ProcessView::ProcessView(QWidget* parent) : QTreeView(parent)
{
}

void ProcessView::updateView()
{
	setAlternatingRowColors(true);
	resizeColumnToContents(0);
	setColumnWidth(1,50);
	setColumnWidth(2,50);
}

void ProcessView::on_Clicked(QModelIndex idx)
{
	qDebug("%d, %d",idx.parent().row(), idx.row());
}

void ProcessView::selectionChanged( const QItemSelection &selected,
                                      const QItemSelection &deselected )
{
	emit itemSelectionChanged(selected, deselected);
	return QTreeView::selectionChanged(selected, deselected);
}

void ProcessView::contextMenuEvent( QContextMenuEvent* event )
{
	Q_UNUSED(event);

	ProcessItem* pcs_item = dynamic_cast<ProcessItem*>(static_cast<ProcessModel*>(this->model())->
	                                          getItem(this->selectionModel()->currentIndex()));
	CondObjectListItem* cond_item =
	        dynamic_cast<CondObjectListItem*>(static_cast<ProcessModel*>(this->model())->
	                                          getItem(this->selectionModel()->currentIndex()));

	if (pcs_item || cond_item)
	{
		QMenu menu;

		if (cond_item)
		{
			QAction* removeCondAction    = menu.addAction("Remove conditions");
			connect(removeCondAction, SIGNAL(triggered()), this, SLOT(removeCondition()));
		}

		if (pcs_item)
		{
			QAction* addCNDAction = menu.addAction("Add FEM Conditions...");
			QAction* removePCSAction    = menu.addAction("Remove process");
			connect(addCNDAction, SIGNAL(triggered()), this, SLOT(addFEMConditions()));
			connect(removePCSAction, SIGNAL(triggered()), this, SLOT(removeProcess()));
		}

		menu.exec(event->globalPos());
	}
}

void ProcessView::removeCondition()
{
	CondObjectListItem* item = dynamic_cast<CondObjectListItem*>(static_cast<ProcessModel*>(this->model())->getItem(this->selectionModel()->currentIndex()));
	
	if (item)
	{
		const FiniteElement::ProcessType pcs_type = static_cast<ProcessItem*>(item->parentItem())->getItem()->getProcessType();
		const FEMCondition::CondType cond_type = item->getType();
		emit conditionsRemoved(pcs_type, cond_type);
	}
}

void ProcessView::removeProcess()
{
	ProcessItem* item = dynamic_cast<ProcessItem*>(static_cast<ProcessModel*>(this->model())->getItem(this->selectionModel()->currentIndex()));
	
	if (item)
	{
		const FiniteElement::ProcessType pcs_type = item->getItem()->getProcessType();
		emit processRemoved(pcs_type);
	}
}

void ProcessView::addFEMConditions()
{
	TreeItem* item = static_cast<ProcessModel*>(model())->getItem(
	        this->selectionModel()->currentIndex());
	emit loadFEMCondFileRequested(item->data(0).toString().toStdString());
}
