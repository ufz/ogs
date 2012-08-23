/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file OGSError.cpp
 *
 * Created on by Karsten Rink
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
void OGSError::box(const QString &e)
{
	box(e, "OpenGeoSys");
}

/**
 * Displays an error in a QMessageBox
 * \param e The error message.
 * \param t The title of the message box
 * \sa QMessageBox
 */
void OGSError::box(const QString &e, const QString &t)
{
	QMessageBox msgBox;
	msgBox.setWindowTitle(t);
	msgBox.setText(e);
	msgBox.exec();
}
