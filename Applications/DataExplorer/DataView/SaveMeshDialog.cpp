/**
 * \file
 * \author Karsten Rink
 * \date   2014-10-27
 * \brief  Implementation of the SaveMeshDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "SaveMeshDialog.h"

#include <QFileDialog>
#include <QSettings>

#include "Base/LastSavedFileDirectory.h"
#include "Base/OGSError.h"
#include "MeshLib/IO/Legacy/MeshIO.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/Mesh.h"

SaveMeshDialog::SaveMeshDialog(MeshLib::Mesh const& mesh, QDialog* parent)
    : QDialog(parent), _mesh(mesh)
{
    setupUi(this);
    this->fileNameEdit->setText(LastSavedFileDirectory::getDir() +
                                QString::fromStdString(_mesh.getName()) +
                                ".vtu");
}

void SaveMeshDialog::on_selectDirButton_clicked()
{
    QString file_type("VTK Unstructured Grid (*.vtu)");
#ifndef NDEBUG
    file_type.append(";;Legacy geometry file (*.msh)");
#endif  // DEBUG
    QSettings settings;
    QString const file_name = QFileDialog::getSaveFileName(
        this,
        "Save mesh as...",
        LastSavedFileDirectory::getDir() +
            QString::fromStdString(_mesh.getName()),
        file_type);

    if (!file_name.isEmpty())
    {
        this->fileNameEdit->setText(file_name);
    }
}

void SaveMeshDialog::on_dataModeBox_currentIndexChanged(int index)
{
    // Disable compression on Ascii
    if (index == 0)
    {
        this->compressionCheckBox->setChecked(false);
        this->compressionCheckBox->setEnabled(false);
        this->compressionLabel->setEnabled(false);
    }
    else
    {
        this->compressionCheckBox->setEnabled(true);
        this->compressionLabel->setEnabled(true);
    }
}

void SaveMeshDialog::accept()
{
    QString const& file_name(this->fileNameEdit->text());
    if (file_name.isEmpty())
    {
        OGSError::box("No file name entered.");
        return;
    }

    QFileInfo fi(file_name);
    if (fi.suffix().toLower() == "vtu")
    {
        int dataMode = this->dataModeBox->currentIndex();
        bool compress(this->compressionCheckBox->isChecked());
        MeshLib::IO::VtuInterface vtkIO(&_mesh, dataMode, compress);
        vtkIO.writeToFile(file_name.toStdString());
    }
    if (fi.suffix().toLower() == "msh")
    {
        MeshLib::IO::Legacy::MeshIO meshIO;
        meshIO.setMesh(&_mesh);
        BaseLib::IO::writeStringToFile(meshIO.writeToString(),
                                       file_name.toStdString());
    }
    LastSavedFileDirectory::setDir(file_name);

    this->done(QDialog::Accepted);
}
