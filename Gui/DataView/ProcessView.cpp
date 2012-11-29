/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file ProcessView.cpp
 *
 * Created on 2010-12-13 by Karsten Rink
 */

#include <QMenu>
#include <QFileDialog>

#include "ProcessItem.h"
#include "CondObjectListItem.h"
#include "CondItem.h"
#include "ProcessModel.h"
#include "ProcessView.h"
#include "FEMConditionSetupDialog.h"
#include "SelectMeshDialog.h"

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

	const QModelIndex idx(this->selectionModel()->currentIndex());
	QMenu menu;

	if (this->isProcessItem(idx))
	{
		QAction* saveCondAction  = menu.addAction("Save FEM Conditions...");
		QAction* removePCSAction = menu.addAction("Remove process");
		connect(saveCondAction, SIGNAL(triggered()), this, SLOT(saveConditions()));
		connect(removePCSAction, SIGNAL(triggered()), this, SLOT(removeProcess()));
	}
	else if (this->isListItem(idx))
	{
		QAction* removeCondAction = menu.addAction("Remove conditions");
		connect(removeCondAction, SIGNAL(triggered()), this, SLOT(removeCondition()));
	}
	else if (this->isConditionItem(idx))
	{
		QAction* editCondAction = menu.addAction("Edit condition");
		connect(editCondAction, SIGNAL(triggered()), this, SLOT(editCondition()));
	}

	menu.exec(event->globalPos());
}

void ProcessView::removeCondition()
{
	CondObjectListItem* item = dynamic_cast<CondObjectListItem*>(static_cast<ProcessModel*>(this->model())->getItem(this->selectionModel()->currentIndex()));

	if (item)
	{
		const FiniteElement::ProcessType pcs_type = static_cast<ProcessItem*>(item->parentItem())->getItem()->getProcessType();
		const FEMCondition::CondType cond_type = item->getType();
		emit conditionsRemoved(pcs_type, "", cond_type);
	}
}

void ProcessView::editCondition()
{
	CondItem* item = dynamic_cast<CondItem*>(static_cast<ProcessModel*>(this->model())->getItem(this->selectionModel()->currentIndex()));

	if (item)
	{
		FEMConditionSetupDialog dlg(*(item->getItem()));
		connect(&dlg, SIGNAL(createFEMCondition(std::vector<FEMCondition*>)), this, SLOT(replaceCondition(std::vector<FEMCondition*>)));
		dlg.exec();
	}
}

void ProcessView::replaceCondition(std::vector<FEMCondition*> conditions)
{
	static_cast<ProcessModel*>(this->model())->replaceCondition(this->selectionModel()->currentIndex(), conditions[0]);
	this->reset();
}

void ProcessView::saveConditions()
{
	emit saveConditionsRequested();
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

bool ProcessView::isProcessItem(const QModelIndex &idx) const
{
	ProcessItem* pcs_item = dynamic_cast<ProcessItem*>(static_cast<ProcessModel*>(this->model())->getItem(idx));
	if (pcs_item) return true;
	return false;
}

bool ProcessView::isListItem(const QModelIndex &idx) const
{
	CondObjectListItem* cond_item = dynamic_cast<CondObjectListItem*>(static_cast<ProcessModel*>(this->model())->getItem(idx));
	if (cond_item) return true;
	return false;
}

bool ProcessView::isConditionItem(const QModelIndex &idx) const
{
	CondObjectListItem* cond_item = dynamic_cast<CondObjectListItem*>(static_cast<ProcessModel*>(this->model())->getItem(idx.parent()));
	if (cond_item) return true;
	return false;
}