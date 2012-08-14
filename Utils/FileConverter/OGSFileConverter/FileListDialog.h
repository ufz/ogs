/**
 * \file FileListDialog.h
 * 2012/04/04 KR Initial implementation
 */

#ifndef FILELISTDIALOG_H
#define FILELISTDIALOG_H

#include "ui_FileList.h"
#include <QDialog>
#include <QStringListModel>

class FileListDialog : public QDialog, private Ui_FileList
{
	Q_OBJECT

public:
	enum FileType {
		GML,	// xml-geometries
		CND,	// xml-fem-conditions
		GLI,	// ascii-geometries
		BC		// ascii-fem-conditions
	};

	FileListDialog(FileType input, FileType output, QWidget* parent = NULL);
	~FileListDialog(void);

private:
	QString getFileTypeString(FileType file_type);

	QStringListModel _allFiles;
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

signals:
	void fileLists(const QStringList, const QString);
};

#endif //FILELISTDIALOG_H
