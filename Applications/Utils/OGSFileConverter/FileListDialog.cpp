/**
 * \file FileListDialog.cpp
 * 2012/04/04 KR Initial implementation
 */

#include "FileListDialog.h"
#include "OGSError.h"

#include "StringTools.h"
#include <QSettings>
#include <QFileInfo>
#include <QLineEdit>

FileListDialog::FileListDialog(FileType input, FileType output, QWidget* parent)
: QDialog(parent), _output_dir(""), _input_file_type(input), _output_file_type(output)
{
	setupUi(this);
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
	QStringList file_names;  

	if (dlg.exec())
	{
		file_names = dlg.selectedFiles();
		if (!file_names.empty()) 
		{
			QStringList list = _allFiles.stringList();
			list.append(file_names);
			_allFiles.setStringList(list);
			QDir dir = QDir(file_names[0]);
			settings.setValue("lastOpenedOgsFileDirectory", dir.absolutePath());
		}
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
	QSettings settings("UFZ", "OpenGeoSys-5");
	QFileInfo fi(settings.value("lastOpenedOgsFileDirectory").toString());
	const QString dirName = QFileDialog::getExistingDirectory(this, "Save to", fi.absolutePath().append("/"));
	this->outputDirEdit->setText(dirName);
}

void FileListDialog::accept()
{
	if (!_allFiles.stringList().empty())
	{
		_output_dir = this->outputDirEdit->text();
		if (!this->outputDirEdit->text().isEmpty() && QDir(_output_dir).exists())
			this->done(QDialog::Accepted);
		else
			OGSError::box("Output directory not found.");
	}
	else
		OGSError::box("No files selected.");
}

void FileListDialog::reject()
{
	this->done(QDialog::Rejected);
}

QString FileListDialog::getFileTypeString(FileType file_type)
{
	if (file_type==GML)		return "OpenGeoSys geometry files (*.gml)";
	else if (file_type==VTU) return "OpenGeoSys mesh files (*.vtu)";
	else if (file_type==GLI) return "GeoSys geometry files (*.gli)";
	else if (file_type==MSH) return "GeoSys mesh files (*.msh)";
	else return "All files (*.*)";
}