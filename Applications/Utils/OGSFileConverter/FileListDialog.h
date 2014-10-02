/**
 * \file   FileListDialog.h
 * \author Karsten Rink
 * \date   2012-04-04
 * \brief  Definition of FileListDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef FILELISTDIALOG_H
#define FILELISTDIALOG_H

#include "ui_FileList.h"
#include <QDialog>
#include <QFileDialog>
#include <QStringListModel>

class FileListDialog : public QDialog, private Ui_FileList
{
	Q_OBJECT

public:
	enum FileType {
		GML,	// xml-geometries
		VTU,	// xml-meshes
		CND,	// xml-fem-conditions
		GLI,	// ascii-geometries
		MSH,	// ascii-meshes
		BC		// ascii-fem-conditions
	};

	FileListDialog(FileType input, FileType output, QWidget* parent = NULL);
	~FileListDialog(void);

	const QStringList getInputFileList() const { return _allFiles.stringList(); };
	const QString getOutputDir() const { return _output_dir; };

private:
	const QString getFileTypeString(FileType file_type) const;

	QStringListModel _allFiles;
	QString _output_dir;
	const FileType _input_file_type;
	const FileType _output_file_type;

private slots:
	void on_addButton_pressed();
	void on_removeButton_pressed();
	void on_browseButton_pressed();

	/// Instructions if the OK-Button has been pressed.
	void accept();

	/// Instructions if the Cancel-Button has been pressed.
	void reject();

};

#endif //FILELISTDIALOG_H
