/**
 * \file DataView.cpp
 * 24/9/2009 LB Initial implementation
 *
 * Implementation of DataView
 */

#include "DataView.h"
#include <QHeaderView>

DataView::DataView( QWidget* parent /*= 0*/ )
: QTreeView(parent)
{
	//verticalHeader()->hide();
	//resizeColumnsToContents();
	//resizeRowsToContents();
}

void DataView::updateView()
{
	setAlternatingRowColors(true);
	size_t nColumns = (this->model() != NULL) ? this->model()->columnCount() : 0;
	for (size_t i=0; i<nColumns; i++)
		resizeColumnToContents(i);
}


/* Selection/Deselection-Functionality -- currently not implemented

QModelIndexList DataView::selectedIndexes() const
{
	return QTreeView::selectedIndexes();
}


void DataView::selectionChanged( const QItemSelection &selected, const QItemSelection &deselected )
{
	emit itemSelectionChanged(selected, deselected);
	return QTreeView::selectionChanged(selected, deselected);
}

void DataView::selectionChangedFromOutside( const QItemSelection &selected, const QItemSelection &deselected )
{
	QItemSelectionModel* selModel = this->selectionModel();

	Q_ASSERT(selModel);

	selModel->blockSignals(true);
	selModel->select(deselected, QItemSelectionModel::Deselect);
	selModel->select(selected, QItemSelectionModel::Select);
	selModel->blockSignals(false);

	Model* model = static_cast<Model*>(this->model());
	//model->setSelectionFromOutside(selected, deselected);

	return QTreeView::selectionChanged(selected, deselected);
}

void DataView::clearSelection()
{
	selectionModel()->clearSelection();
}
*/
