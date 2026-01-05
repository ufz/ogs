// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <QDir>
#include <QFileInfo>
#include <QSettings>
#include <QString>

class LastSavedFileDirectory
{
public:
    /// Returns the directory last used for saving a file
    static const QString getDir()
    {
        QSettings settings;
        return settings.value("lastSavedFileDirectory").toString();
    }

    /// Sets the directory last used for saving a file
    static void setDir(const QString &path)
    {
        QFileInfo fi(path);
        QDir dir = QDir(fi.absolutePath());
        QSettings settings;
        settings.setValue("lastSavedFileDirectory", dir.absolutePath() + "/");
    }

};
