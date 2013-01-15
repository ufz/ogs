/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the OGSError class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
