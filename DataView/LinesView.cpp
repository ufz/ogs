/**
 * \file QOgsDataView.cpp
 * 24/9/2009 LB Initial implementation
 * 
 * Implementation of QOgsDataView
 */

#include "QOgsLinesView.h"

#include <QHeaderView>

QOgsDataView::QOgsDataView( QWidget* parent /*= 0*/ )
: QTableView(parent)
{
	verticalHeader()->hide();
	resizeColumnsToContents();
	resizeRowsToContents();
}

void QOgsDataView::selectionChanged( const QItemSelection &selected, const QItemSelection &deselected )
{
	emit itemSelectionChanged(selected, deselected);
	return QTableView::selectionChanged(selected, deselected);
}

void QOgsDataView::selectionChangedFromOutside( const QItemSelection &selected, const QItemSelection &deselected )
{
	QAbstractItemModel* model = this->model();
	QItemSelection deselection(model->index(0, 0), model->index(model->rowCount()-1, 0));
	QItemSelectionModel* selModel = this->selectionModel();

	selModel->blockSignals(true);
	selModel->select(deselection, QItemSelectionModel::Deselect);
	selModel->select(selected, QItemSelectionModel::Select);
	selModel->blockSignals(false);

	return QTableView::selectionChanged(selected, deselected);
}

QModelIndexList QOgsDataView::selectedIndexes() const
{
	return QTableView::selectedIndexes();
}