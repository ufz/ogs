/**
 * \file OGSError.cpp
 * KR Initial implementation
 */
#include "OGSError.h"

#include <QMessageBox>
#include <QString>

OGSError::OGSError()
{}

OGSError::~OGSError()
{}

/**
 * Displays an error in a QMessageBox
 * \param e The error message.
 */
void OGSError::box(QString e)
{
	box(e, "OpenGeoSys");
}

/**
 * Displays an error in a QMessageBox
 * \param e The error message.
 * \param t The title of the message box
 * \sa QMessageBox
 */
void OGSError::box(QString e, QString t)
{
	QMessageBox msgBox;
	msgBox.setWindowTitle(t);
	msgBox.setText(e);
	msgBox.exec();
}
