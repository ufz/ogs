/**
 * \file QueryResultsDialog.cpp
 * KR Initial implementation
 */

#include "OGSError.h"
#include "QueryResultsDialog.h"

/**
 * Constructor.
 */
QueryResultsDialog::QueryResultsDialog(QDialog* parent) : QDialog(parent)
{
	setupUi(this);
}

QueryResultsDialog::~QueryResultsDialog()
{
}

void QueryResultsDialog::setView(QSqlQueryModel* model)
{
	queryView->setModel(model);
	setViewProperties();
}

void QueryResultsDialog::setViewProperties()
{
	queryView->setColumnHidden(0, true); //hide ID column
	queryView->setSelectionMode(QAbstractItemView::SingleSelection);
	queryView->setSelectionBehavior(QAbstractItemView::SelectRows);
	queryView->setEditTriggers(QAbstractItemView::NoEditTriggers);
	queryView->resizeColumnsToContents();

	QHeaderView* header = queryView->horizontalHeader();
	header->setStretchLastSection(true);
}

void QueryResultsDialog::on_openButton_clicked()
{
	QItemSelectionModel* selection = queryView->selectionModel();

	if (selection->hasSelection())
	{
		queryView->selectColumn(0);
		QModelIndexList indexes = selection->selectedIndexes();
		int listID = (queryView->model()->data(indexes.first())).toInt();

		emit listSelected(listID);

		this->done(QDialog::Accepted);
	}
	else
		OGSError::box("No data selected.");
}

void QueryResultsDialog::on_cancelButton_clicked()
{
	this->done(QDialog::Rejected);
}
