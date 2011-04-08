/**
 * \file ConditionView.cpp
 * 2010/12/13 KR Initial implementation
 */

#include <QMenu>

#include "ConditionView.h"
#include "ConditionModel.h"
#include "CondObjectListItem.h"
#include "CondItem.h"


ConditionView::ConditionView(QWidget* parent) : QTreeView(parent)
{
}

void ConditionView::updateView()
{
	setAlternatingRowColors(true);
	resizeColumnToContents(0);
	setColumnWidth(1,50);
	setColumnWidth(2,50);
}

void ConditionView::on_Clicked(QModelIndex idx)
{
	qDebug("%d, %d",idx.parent().row(), idx.row());
}

void ConditionView::selectionChanged( const QItemSelection &selected, const QItemSelection &deselected )
{
	emit itemSelectionChanged(selected, deselected);
	return QTreeView::selectionChanged(selected, deselected);
}

void ConditionView::contextMenuEvent( QContextMenuEvent* event )
{
	Q_UNUSED(event);
	
	CondObjectListItem* item = dynamic_cast<CondObjectListItem*>(static_cast<ConditionModel*>(this->model())->getItem(this->selectionModel()->currentIndex()));
	if (item)
	{
		QMenu menu;
		QAction* removeAction    = menu.addAction("Remove");
		connect(removeAction, SIGNAL(triggered()), this, SLOT(removeCondition()));
		menu.exec(event->globalPos());
	}

}

void ConditionView::removeCondition()
{
	CondObjectListItem* item = dynamic_cast<CondObjectListItem*>(static_cast<ConditionModel*>(this->model())->getItem(this->selectionModel()->currentIndex()));
	QString geo_name = item->parentItem()->data(0).toString();
	FEMCondition::CondType type = item->getType();
	emit conditionsRemoved(geo_name, type);
}
/*
void ConditionView::removeAllConditions()
{
	ConditionModel* model = static_cast<ConditionModel*>(this->model());

	for (size_t j=0; j<3; j++)
	{
		QModelIndex parentIndex = model->index(j, 0, QModelIndex());
		int nChildren = model->getItem(parentIndex)->childCount();
		for (int i=nChildren; i>=0; i--)
			emit requestConditionRemoval(model->index(i, 0, parentIndex));
	}
}
*/
