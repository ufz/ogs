/**
 * \file   FileListDialog.h
 * \author Karsten Rink
 * \date   2012-04-04
 * \brief  Definition of FileListDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ui_FileList.h"
#include <QDialog>
#include <QFileDialog>
#include <QStringListModel>

enum class FileType {
    GML,    // xml-geometries
    VTU,    // xml-meshes
    GLI,    // ascii-geometries
    MSH,    // ascii-meshes
};

/**
 * A list of selected files selected for conversion incl. an output directory.
 */
class FileListDialog : public QDialog, private Ui_FileList
{
    Q_OBJECT

public:
    /// Constructor
    FileListDialog(FileType input, FileType output, QWidget* parent = nullptr);
    /// Destructor
    ~FileListDialog(void) override;

    /// Returns list of all selected files
    const QStringList getInputFileList() const { return _allFiles.stringList(); };
    /// Returns selected output directory
    const QString getOutputDir() const { return _output_dir; };

private:
    /// Returns a string for the given file type enum
    const QString getFileTypeString(FileType file_type) const;
    /// Display a warning for vtu- to msh-conversion
    void displayWarningLabel() const;

    QStringListModel _allFiles;
    QString _output_dir;
    const FileType _input_file_type;
    const FileType _output_file_type;

private slots:
    void on_addButton_pressed();
    void on_removeButton_pressed();
    void on_browseButton_pressed();

    /// Instructions if the OK-Button has been pressed.
    void accept() override;

    /// Instructions if the Cancel-Button has been pressed.
    void reject() override;
};
