/**
 * \file QueryResultsDialog.h
 * KR Initial implementation
 */

#ifndef QUERYRESULTSDIALOG_H
#define QUERYRESULTSDIALOG_H

#include "ui_DatabaseResultView.h"
#include <QSqlQueryModel>
#include <QDialog>

/**
 * \brief A Dialog for displaying the results of a database query in a table.
 */
class QueryResultsDialog : public QDialog, private Ui_DatabaseResultsView
{
	Q_OBJECT

public:
	QueryResultsDialog(QDialog* parent = 0);
	~QueryResultsDialog(void);

	/// Sets up the view.
	void setView(QSqlQueryModel* model);

private:
	/// Sets the properties of the view.
	void setViewProperties();

private slots:
	/// Instructions if the Cancel-button is clicked.
	void on_cancelButton_clicked();

	/// Instructions if the Open-button is clicked.
	void on_openButton_clicked();

signals:
	void listSelected(int listID);
};

#endif //QUERYRESULTSDIALOG_H
