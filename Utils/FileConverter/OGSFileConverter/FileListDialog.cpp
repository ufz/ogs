/**
 * \file FileListDialog.cpp
 * 2012/04/04 KR Initial implementation
 */

#include "FileListDialog.h"

#include "StringTools.h"
#include <QFileDialog>
#include <QSettings>
#include <QFileInfo>

FileListDialog::FileListDialog(FileType input, FileType output, QWidget* parent)
: QDialog(parent), _input_file_type(input), _output_file_type(output)
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
	QString fileName = QFileDialog::getOpenFileName( this, "Select data file to open",
		                                             settings.value("lastOpenedOgsFileDirectory").toString(), 
													 this->getFileTypeString(_input_file_type));
	
	if (!fileName.isEmpty())
	{
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedOgsFileDirectory", dir.absolutePath());
		QStringList list = _allFiles.stringList();
		list.append(fileName);
		_allFiles.setStringList(list);
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
	QString guess_name("");
	if (!_allFiles.stringList().empty())
		guess_name = QString::fromStdString(BaseLib::getFileNameFromPath(_allFiles.stringList().at(0).toStdString()));
	QSettings settings("UFZ", "OpenGeoSys-5");
	QFileInfo fi(settings.value("lastOpenedOgsFileDirectory").toString());
	QString fileName = QFileDialog::getSaveFileName( this, "Save as",
													 fi.absolutePath().append("/").append(guess_name), 
													 this->getFileTypeString(_output_file_type));
	
	if (!fileName.isEmpty())
	{
		QDir dir = QDir(fileName);
		settings.setValue("lastOpenedOgsFileDirectory", dir.absolutePath());
		this->lineEdit->setText(fileName);
	}
}

void FileListDialog::accept()
{
	emit fileLists(_allFiles.stringList(), lineEdit->text());
	this->done(QDialog::Accepted);
}

void FileListDialog::reject()
{
	this->done(QDialog::Rejected);
}

QString FileListDialog::getFileTypeString(FileType file_type)
{
	if (file_type==GML)		return "OpenGeoSys geometry files (*.gml)";
	else if (file_type==CND) return "OpenGeoSys condition files (*.cnd)";
	else if (file_type==GLI) return "GeoSys geometry files (*.gli)";
	else if (file_type==BC)	return "GeoSys boundary condition files (*.bc);; GeoSys initial condition files (*.ic);; GeoSys source term files (*.st)";
	else return "All files (*.*)";
}