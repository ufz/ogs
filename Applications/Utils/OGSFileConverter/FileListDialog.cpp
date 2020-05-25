/**
 * \file
 * \author Karsten Rink
 * \date   2012-04-04
 * \brief  Implementation of FileListDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "FileListDialog.h"

#include <QSettings>
#include <QFileInfo>
#include <QFont>
#include <QLineEdit>

#include "Applications/DataExplorer/Base/OGSError.h"
#include "BaseLib/StringTools.h"

FileListDialog::FileListDialog(FileType input, FileType output, QWidget* parent)
: QDialog(parent), output_dir_(""), input_file_type_(input), output_file_type_(output)
{
    setupUi(this);

    if (input_file_type_ == FileType::VTU && output_file_type_ == FileType::MSH)
    {
        displayWarningLabel();
    }

    this->listView->setModel(&allFiles_);
}

FileListDialog::~FileListDialog() = default;

void FileListDialog::on_addButton_pressed()
{
    QSettings settings("UFZ", "OpenGeoSys-5");
    QFileDialog dlg(this);
    dlg.setDirectory(settings.value("lastOpenedOgsFileDirectory").toString());
    dlg.setFileMode(QFileDialog::ExistingFiles);
    dlg.setNameFilter(this->getFileTypeString(input_file_type_));

    if (dlg.exec())
    {
        QStringList const file_names = dlg.selectedFiles();

        if (file_names.empty())
        {
            return;
        }

        QStringList list = allFiles_.stringList();
        list.append(file_names);
        allFiles_.setStringList(list);
        QDir const dir = QDir(file_names[0]);
        settings.setValue("lastOpenedOgsFileDirectory", dir.absolutePath());
    }
}

void FileListDialog::on_removeButton_pressed()
{
    QModelIndexList selected = this->listView->selectionModel()->selectedIndexes();
    for (auto& item : selected)
    {
        this->allFiles_.removeRow(item.row());
    }
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
    if (allFiles_.stringList().empty())
    {
        OGSError::box("No files selected.");
        return;
    }

    output_dir_ = this->outputDirEdit->text();
    if (!this->outputDirEdit->text().isEmpty() && QDir(output_dir_).exists())
    {
        this->done(QDialog::Accepted);
    }
    else
    {
        OGSError::box("Output directory not found.");
    }
}

void FileListDialog::reject()
{
    this->done(QDialog::Rejected);
}

QString FileListDialog::getFileTypeString(FileType file_type) const
{
    if (file_type == FileType::GML)
    {
        return "OpenGeoSys geometry files (*.gml)";
    }
    if (file_type == FileType::VTU)
    {
        return "OpenGeoSys mesh files (*.vtu)";
    }
    if (file_type == FileType::GLI)
    {
        return "GeoSys geometry files (*.gli)";
    }
    if (file_type == FileType::MSH)
    {
        return "GeoSys mesh files (*.msh)";
    }
    return "All files (*.*)";
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
