/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the OGSError class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#include "Base/OGSError.h"

#include <QMessageBox>
#include <QString>

OGSError::OGSError() = default;

OGSError::~OGSError() = default;

void OGSError::box(const QString& e)
{
    box(e, "OpenGeoSys");
}

void OGSError::box(const QString& e, const QString& t)
{
    QMessageBox msgBox;
    msgBox.setWindowTitle(t);
    msgBox.setText(e);
    msgBox.exec();
}

bool OGSError::question(const QString& e, const QString& t)
{
    QMessageBox msgBox;
    msgBox.setWindowTitle(t);
    msgBox.setText(e);
    msgBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
    msgBox.setDefaultButton(QMessageBox::Cancel);

    return msgBox.exec() == QMessageBox::Ok;
}
