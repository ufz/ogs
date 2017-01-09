/**
 * \file   FileListDialog.cpp
 * \author Karsten Rink
 * \date   2012-04-04
 * \brief  Implementation of FileListDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FileListDialog.h"
#include "OGSError.h"

#include "StringTools.h"
#include <QSettings>
#include <QFileInfo>
#include <QFont>
#include <QLineEdit>

FileListDialog::FileListDialog(FileType input, FileType output, QWidget* parent)
: QDialog(parent), _output_dir(""), _input_file_type(input), _output_file_type(output)
{
    setupUi(this);

    if (_input_file_type == FileType::VTU && _output_file_type == FileType::MSH)
        displayWarningLabel();

    this->listView->setModel(&_allFiles);
}

FileListDialog::~FileListDialog()
{
}

void FileListDialog::on_addButton_pressed()
{
    QSettings settings("UFZ", "OpenGeoSys-5");
    QFileDialog dlg(this);
    dlg.setDirectory(settings.value("lastOpenedOgsFileDirectory").toString());
    dlg.setFileMode(QFileDialog::ExistingFiles);
    dlg.setNameFilter(this->getFileTypeString(_input_file_type));

    if (dlg.exec())
    {
        QStringList const file_names = dlg.selectedFiles();

        if (file_names.empty())
            return;

        QStringList list = _allFiles.stringList();
        list.append(file_names);
        _allFiles.setStringList(list);
        QDir const dir = QDir(file_names[0]);
        settings.setValue("lastOpenedOgsFileDirectory", dir.absolutePath());
    }
}

void FileListDialog::on_removeButton_pressed()
{
    QModelIndexList selected = this->listView->selectionModel()->selectedIndexes();
    for (QModelIndexList::iterator it = selected.begin(); it != selected.end(); ++it)
        this->_allFiles.removeRow(it->row());
}

void FileListDialog::on_browseButton_pressed()
{
    QSettings const settings("UFZ", "OpenGeoSys-5");
    QFileInfo const fi(settings.value("lastOpenedOgsFileDirectory").toString());
    QString const dirName = QFileDialog::getExistingDirectory(this, "Save to", fi.absolutePath().append("/"));
    this->outputDirEdit->setText(dirName);
}

void FileListDialog::accept()
{
    if (_allFiles.stringList().empty())
    {
        OGSError::box("No files selected.");
        return;
    }

    _output_dir = this->outputDirEdit->text();
    if (!this->outputDirEdit->text().isEmpty() && QDir(_output_dir).exists())
        this->done(QDialog::Accepted);
    else
        OGSError::box("Output directory not found.");
}

void FileListDialog::reject()
{
    this->done(QDialog::Rejected);
}

const QString FileListDialog::getFileTypeString(FileType file_type) const
{
    if      (file_type==FileType::GML) return "OpenGeoSys geometry files (*.gml)";
    else if (file_type==FileType::VTU) return "OpenGeoSys mesh files (*.vtu)";
    else if (file_type==FileType::GLI) return "GeoSys geometry files (*.gli)";
    else if (file_type==FileType::MSH) return "GeoSys mesh files (*.msh)";
    else return "All files (*.*)";
}

void FileListDialog::displayWarningLabel() const
{
    this->warningLabel->setFixedWidth(300);
    this->warningLabel->setFixedHeight(40);
    QFont font = this->warningLabel->font();
    font.setPointSize(9);
    font.setBold(true);
    warningLabel->setFont(font);
    this->warningLabel->setStyleSheet("QLabel { color : red; }");
    this->warningLabel->setText("Warning: all scalar information except<br />MaterialIDs will be lost!");
}

