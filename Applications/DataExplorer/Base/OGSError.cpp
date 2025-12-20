// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
