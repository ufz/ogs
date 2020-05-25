/**
 * \file
 * \author Lars Bilke
 * \date   2009-11-05
 * \brief  Implementation of the RecentFiles class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// ** INCLUDES **
#include "RecentFiles.h"

#include <QFileInfo>
#include <QSettings>
#include <utility>

RecentFiles::RecentFiles(QObject* parent, const char* slot,
                         QString settingsName)
    : QObject(parent), settingsName_(std::move(settingsName))
{
    filesMenu_ = new QMenu(tr("Recent files"));
    for (auto& fileAction : fileActions_)
    {
        fileAction = new QAction(this);
        fileAction->setVisible(false);
        connect(fileAction, SIGNAL(triggered()), parent, slot);
        filesMenu_->addAction(fileAction);
    }
    updateRecentFileActions();
}

RecentFiles::~RecentFiles()
{
    delete filesMenu_;
}

QMenu* RecentFiles::menu()
{
    return filesMenu_;
}
void RecentFiles::setCurrentFile( const QString& filename )
{
    currentFile_ = filename;

    QSettings settings;
    QStringList files = settings.value(settingsName_).toStringList();
    files.removeAll(filename);
    files.prepend(filename);
    while (files.size() > maxFiles_)
    {
        files.removeLast();
    }

    settings.setValue("recentFileList", files);

    updateRecentFileActions();
}
void RecentFiles::updateRecentFileActions()
{
    QSettings settings;
    QStringList files = settings.value(settingsName_).toStringList();

    int numFiles = qMin(files.size(), static_cast<int>(maxFiles_));

    for (int i = 0; i < numFiles; ++i)
    {
        QString text = tr("&%1 %2").arg(i + 1).arg(strippedName(files[i]));
        fileActions_[i]->setText(text);
        fileActions_[i]->setData(files[i]);
        fileActions_[i]->setVisible(true);
    }

    for (int i = numFiles; i < maxFiles_; ++i)
    {
        fileActions_[i]->setVisible(false);
    }
}

QString RecentFiles::strippedName( const QString& fullFileName )
{
    return QFileInfo(fullFileName).fileName();
}
