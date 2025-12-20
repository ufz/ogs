// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

// ** INCLUDES **
#include "RecentFiles.h"

#include <QFileInfo>
#include <QSettings>
#include <utility>

RecentFiles::RecentFiles(QObject* parent, const char* slot,
                         QString settingsName)
    : QObject(parent), _settingsName(std::move(settingsName))
{
    _filesMenu = new QMenu(tr("Recent files"));
    for (auto& fileAction : _fileActions)
    {
        fileAction = new QAction(this);
        fileAction->setVisible(false);
        connect(fileAction, SIGNAL(triggered()), parent, slot);
        _filesMenu->addAction(fileAction);
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
void RecentFiles::setCurrentFile(const QString& filename)
{
    _currentFile = filename;

    QSettings settings;
    QStringList files = settings.value(_settingsName).toStringList();
    files.removeAll(filename);
    files.prepend(filename);
    while (files.size() > _maxFiles)
    {
        files.removeLast();
    }

    settings.setValue("recentFileList", files);

    updateRecentFileActions();
}
void RecentFiles::updateRecentFileActions()
{
    QSettings settings;
    QStringList files = settings.value(_settingsName).toStringList();

    int numFiles = qMin(files.size(), static_cast<int>(_maxFiles));

    for (int i = 0; i < numFiles; ++i)
    {
        QString text = tr("&%1 %2").arg(i + 1).arg(strippedName(files[i]));
        _fileActions[i]->setText(text);
        _fileActions[i]->setData(files[i]);
        _fileActions[i]->setVisible(true);
    }

    for (int i = numFiles; i < _maxFiles; ++i)
    {
        _fileActions[i]->setVisible(false);
    }
}

QString RecentFiles::strippedName(const QString& fullFileName)
{
    return QFileInfo(fullFileName).fileName();
}
