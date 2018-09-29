/**
 * \file   LastSavedFileDirectory.h
 * \author Karsten Rink
 * \date   2014-02-04
 * \brief  Manages the last directory used for saving a file
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
