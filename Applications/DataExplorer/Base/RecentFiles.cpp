/**
 * \file
 * \author Lars Bilke
 * \date   2009-11-05
 * \brief  Implementation of the RecentFiles class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "RecentFiles.h"

#include <QFileInfo>
#include <QSettings>

RecentFiles::RecentFiles(  QObject* parent, const char* slot, QString settingsName)
    : QObject(parent), _settingsName(settingsName)
{
    _filesMenu = new QMenu(tr("Recent files"));
    for (int i = 0; i < _maxFiles; i++)
    {
        _fileActions[i] = new QAction(this);
        _fileActions[i]->setVisible(false);
        connect(_fileActions[i], SIGNAL(triggered()), parent, slot);
        _filesMenu->addAction(_fileActions[i]);
    }
    updateRecentFileActions();
}

RecentFiles::~RecentFiles()
{
    delete _filesMenu;
}

QMenu* RecentFiles::menu()
{
    return _filesMenu;
}
void RecentFiles::setCurrentFile( const QString& filename )
{
    _currentFile = filename;

    QSettings settings;
    QStringList files = settings.value(_settingsName).toStringList();
    files.removeAll(filename);
    files.prepend(filename);
    while (files.size() > _maxFiles)
        files.removeLast();

    settings.setValue("recentFileList", files);

    updateRecentFileActions();
}
void RecentFiles::updateRecentFileActions()
{
    QSettings settings;
    QStringList files = settings.value(_settingsName).toStringList();

    int numFiles = qMin(files.size(), (int)_maxFiles);

    for (int i = 0; i < numFiles; ++i)
    {
        QString text = tr("&%1 %2").arg(i + 1).arg(strippedName(files[i]));
        _fileActions[i]->setText(text);
        _fileActions[i]->setData(files[i]);
        _fileActions[i]->setVisible(true);
    }

    for (int i = numFiles; i < _maxFiles; ++i)
        _fileActions[i]->setVisible(false);
}

QString RecentFiles::strippedName( const QString& fullFileName )
{
    return QFileInfo(fullFileName).fileName();
}
