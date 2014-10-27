/**
 * \file   SaveMeshDialog.cpp
 * \author Karsten Rink
 * \date   2014-10-27
 * \brief  Implementation of the SaveMeshDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SaveMeshDialog.h"

#include <QFileDialog>
#include <QSettings>

// ThirdParty/logog
#include "logog/include/logog.hpp"

#include "FileIO/VtkIO/VtuInterface.h"
#include "Legacy/MeshIO.h"

#include "Mesh.h"
#include "OGSError.h"
#include "LastSavedFileDirectory.h"

SaveMeshDialog::SaveMeshDialog(MeshLib::Mesh const& mesh, QDialog* parent)
	: _mesh(mesh), QDialog(parent)
{
	setupUi(this);
	this->fileNameEdit->setText(LastSavedFileDirectory::getDir() + QString::fromStdString(_mesh.getName()) + ".vtu");

	// set default data mode to appended because of ascii bug (see VtuInterface documentation)
	this->dataModeBox->setCurrentIndex(2);
}

void SaveMeshDialog::on_selectDirButton_clicked()
{
	QString file_type ("VTK Unstructured Grid (*.vtu)");
#ifndef NDEBUG
	file_type.append(";;Legacy geometry file (*.msh)");
#endif // DEBUG
	QSettings settings;
	QString const file_name = QFileDialog::getSaveFileName(this, 
		"Save mesh as...",
		LastSavedFileDirectory::getDir() + QString::fromStdString(_mesh.getName()), 
		file_type);

	if (!file_name.isEmpty())
		this->fileNameEdit->setText(file_name);
}

void SaveMeshDialog::accept()
{
	QString const& file_name (this->fileNameEdit->text());
	if (file_name.isEmpty())
	{
		OGSError::box("No file name entered.");
		return;
	}

	QFileInfo fi(file_name);
	if (fi.suffix().toLower() == "vtu")
	{
		bool append (this->dataModeBox->currentIndex() == 2);
		bool binary (this->dataModeBox->currentIndex() > 0);
		bool compress (this->compressionBox->currentIndex() > 0);
		FileIO::VtuInterface vtkIO(&_mesh, binary, append, compress);
		vtkIO.writeToFile(file_name.toStdString().c_str());
	}
	if (fi.suffix().toLower() == "msh")
	{
		FileIO::Legacy::MeshIO meshIO;
		meshIO.setMesh(&_mesh);
		meshIO.writeToFile(file_name.toStdString().c_str());
	}
	LastSavedFileDirectory::setDir(file_name);
			
	this->done(QDialog::Accepted);	
}

