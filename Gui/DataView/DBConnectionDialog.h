/**
 * \file DBConnectionDialog.h
 * KR Initial implementation
 */

#ifndef DBCONNECTIONDIALOG_H
#define DBCONNECTIONDIALOG_H

#include "ui_DBConnection.h"
#include <QSqlQueryModel>
#include <QDialog>

/**
 * \brief A dialog window for settung up a database connection
 */
class DBConnectionDialog : public QDialog, private Ui_DBConnectionDialog
{
	Q_OBJECT

public:
	DBConnectionDialog(QDialog* parent = 0);
	~DBConnectionDialog(void);

private slots:
	void accept();
	void reject();

signals:
	void connectionRequested(QString, QString, QString, QString, QString);
};

#endif //DBCONNECTIONDIALOG_H
