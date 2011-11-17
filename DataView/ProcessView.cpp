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

	CondObjectListItem* item =
	        dynamic_cast<CondObjectListItem*>(static_cast<ProcessModel*>(this->model())->
	                                          getItem(this->selectionModel()->currentIndex()));
	if (item)
	{
		QMenu menu;
		QAction* removeAction    = menu.addAction("Remove");
		connect(removeAction, SIGNAL(triggered()), this, SLOT(removeCondition()));
		menu.exec(event->globalPos());
	}
}

void ProcessView::removeCondition()
{
	CondObjectListItem* item =
	        dynamic_cast<CondObjectListItem*>(static_cast<ProcessModel*>(this->model())->
	                                          getItem(this->selectionModel()->currentIndex()));
	QString process_name = item->parentItem()->data(0).toString();
	FEMCondition::CondType type = item->getType();
	emit conditionsRemoved(process_name, type);
}
/*
   void ProcessView::removeAllConditions()
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
