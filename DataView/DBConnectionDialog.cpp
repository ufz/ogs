/**
 * \file DBConnectionDialog.cpp
 * KR Initial implementation
 */

#include <QSettings>
#include "DBConnectionDialog.h"


/// Constructor
DBConnectionDialog::DBConnectionDialog(QDialog* parent) : QDialog(parent)
{
	setupUi(this);

	int idx=0;
	QSettings settings("UFZ", "OpenGeoSys-5");

	if (!settings.value("DBProtocol", "").toString().isEmpty())
	{
		for (int i=0; i<driverBox->count(); i++)
		{
			if (driverBox->itemText(i).startsWith(settings.value("DBProtocol", "").toString()))
				idx = i;
		}
	}

	driverBox->setCurrentIndex(idx);
	hostnameLine->setText(settings.value("DBHost",     "").toString());
	dbnameLine->setText(settings.value("DBName",       "").toString());
	usernameLine->setText(settings.value("DBUser",     "").toString());
	passwordLine->setText(settings.value("DBPass",     "").toString());
}

DBConnectionDialog::~DBConnectionDialog()
{
}

/// Instructions if the OK-Button has been pressed.
void DBConnectionDialog::accept()
{
	QString protocol = (driverBox->currentText()).left((driverBox->currentText()).indexOf(" "));
	QString hostname = hostnameLine->text();
	QString dbname   = dbnameLine->text();
	QString user     = usernameLine->text();
	QString pass     = passwordLine->text();

	emit connectionRequested(protocol, hostname, dbname, user, pass);
	this->done(QDialog::Accepted);
}

/// Instructions if the Cancel-Button has been pressed.
void DBConnectionDialog::reject()
{
	this->done(QDialog::Rejected);
}
