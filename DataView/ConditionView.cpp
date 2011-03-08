/**
 * \file ConditionView.cpp
 * 2010/12/13 KR Initial implementation
 */

#include <QMenu>

#include "ConditionView.h"
#include "ConditionModel.h"
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
//	QModelIndex index = this->selectionModel()->currentIndex();
//	CondItem* item = static_cast<CondItem*>(index.internalPointer());

}

void ConditionView::removeCondition()
{
	TreeItem* item = static_cast<ConditionModel*>(model())->getItem(this->selectionModel()->currentIndex());
	emit conditionRemoved((item->data(0).toString()).toStdString());
}

