/**
 * \file
 * \date   2023-05-11
 * \brief  Implementation of the Vtu2GridDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Vtu2GridDialog.h"

#include <vtkXMLUnstructuredGridWriter.h>

#include <QStringList>
#include <QStringListModel>
#include <string>

#include "Base/StrictDoubleValidator.h"
#include "GeoLib/AABB.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Vtk/VtkMappedMeshSource.h"
#include "MeshModel.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"
#include "MeshToolsLib/MeshGenerators/VoxelGridFromMesh.h"

Vtu2GridDialog::Vtu2GridDialog(MeshModel& mesh_model, QDialog* parent)
    : QDialog(parent), _mesh_model(mesh_model)
{
    setupUi(this);
    QStringList MeshList;

    for (int model_index = 0; model_index < mesh_model.rowCount();
         ++model_index)
    {
        auto const* mesh = mesh_model.getMesh(mesh_model.index(model_index, 0));
        MeshList.append(QString::fromStdString(mesh->getName()));
    }

    if (MeshList.empty())
    {
        MeshList.append("[No Mesh available.]");
        this->expectedVoxelLabel->setText("Expected Voxel: undefined");
    }

    _allMeshes.setStringList(MeshList);
    this->meshListBox->addItems(_allMeshes.stringList());
    this->xlineEdit->setFocus();
}

void Vtu2GridDialog::updateExpectedVoxel()
{
    QString const xin = this->xlineEdit->text();
    QString const yin = this->ylineEdit->text();
    QString const zin = this->zlineEdit->text();
    bool ok;
    double const xinput = xin.toDouble();
    double const yinput = (yin.toDouble(&ok)) ? yin.toDouble() : xin.toDouble();
    double const zinput = (zin.toDouble(&ok)) ? zin.toDouble() : xin.toDouble();

    if (_allMeshes.stringList()[0] == "[No Mesh available.]")
    {
        this->expectedVoxelLabel->setText("approximated Voxel: undefined");
        return;
    }
    if (xin.isEmpty() || xinput == 0)
    {
        this->expectedVoxelLabel->setText("approximated Voxel: undefined");
        return;
    }

    auto* const _mesh(
        _mesh_model.getMesh(this->meshListBox->currentText().toStdString()));
    auto const& nodes = _mesh->getNodes();
    GeoLib::AABB const aabb(nodes.cbegin(), nodes.cend());

    auto const min = aabb.getMinPoint();
    auto const max = aabb.getMaxPoint();

    double const expectedVoxel = (max[0] - min[0]) * (max[1] - min[1]) *
                                 (max[2] - min[2]) / xinput / yinput / zinput;

    int const exponent = std::floor(std::log10(abs(expectedVoxel)));
    this->expectedVoxelLabel->setText(
        "approximated Voxel = " +
        QString::number(std::round(expectedVoxel / std::pow(10, exponent))) +
        " x 10^" + QString::number(exponent));
}

void Vtu2GridDialog::on_xlineEdit_textChanged()
{
    updateExpectedVoxel();
}

void Vtu2GridDialog::on_ylineEdit_textChanged()
{
    updateExpectedVoxel();
}

void Vtu2GridDialog::on_zlineEdit_textChanged()
{
    updateExpectedVoxel();
}

void Vtu2GridDialog::accept()
{
    using namespace MeshToolsLib::MeshGenerator;
    if (this->meshListBox->currentText().toStdString() ==
        "[No Mesh available.]")
    {
        OGSError::box(
            "Please specify the input meshes. It has to be a 3D mesh.");
        return;
    }

    QString const xin = this->xlineEdit->text();
    QString const yin = this->ylineEdit->text();
    QString const zin = this->zlineEdit->text();

    bool ok;
    if (!xin.toDouble(&ok))
    {
        OGSError::box(
            "At least the x-length of a voxel must be specified.\n If "
            "y-/z-input "
            "are not specified, equal to 0, or not a real number, they are "
            "treated as "
            "the x-input.");
        return;
    }
    double const xinput = xin.toDouble();
    double const yinput = (yin.toDouble(&ok)) ? yin.toDouble() : xin.toDouble();
    double const zinput = (zin.toDouble(&ok)) ? zin.toDouble() : xin.toDouble();
    std::array<double, 3> const cellsize = {xinput, yinput, zinput};

    auto _mesh(
        _mesh_model.getMesh(this->meshListBox->currentText().toStdString()));

    if (_mesh->MeshLib::Mesh::getDimension() < 3)
    {
        OGSError::box("The dimension of the mesh has to be 3.");
        return;
    }

    vtkNew<MeshLib::VtkMappedMeshSource> vtkSource;
    vtkSource->SetMesh(_mesh);
    vtkSource->Update();
    vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSource->GetOutput();

    double* const bounds = mesh->GetBounds();
    MathLib::Point3d const min(
        std::array<double, 3>{bounds[0], bounds[2], bounds[4]});
    MathLib::Point3d const max(
        std::array<double, 3>{bounds[1], bounds[3], bounds[5]});
    std::array<double, 3> ranges = {max[0] - min[0], max[1] - min[1],
                                    max[2] - min[2]};
    if (ranges[0] < 0 || ranges[1] < 0 || ranges[2] < 0)
    {
        OGSError::box(
            "The range (max-min of the bounding box) is not allowed to be < 0");
    }
    std::array<std::size_t, 3> const dims =
        VoxelGridFromMesh::getNumberOfVoxelPerDimension(ranges, cellsize);
    std::unique_ptr<MeshLib::Mesh> grid(
        generateRegularHexMesh(dims[0], dims[1], dims[2], cellsize[0],
                               cellsize[1], cellsize[2], min, "grid"));

    std::vector<int> const tmp_ids =
        VoxelGridFromMesh::assignCellIds(mesh, min, dims, cellsize);
    std::vector<int>& cell_ids =
        *grid->getProperties().createNewPropertyVector<int>(
            VoxelGridFromMesh::cell_id_name, MeshLib::MeshItemType::Cell, 1);
    std::copy(tmp_ids.cbegin(), tmp_ids.cend(), std::back_inserter(cell_ids));

    if (!VoxelGridFromMesh::removeUnusedGridCells(mesh, grid))
    {
        return;
    }

    VoxelGridFromMesh::mapMeshArraysOntoGrid(mesh, grid);

    if (mesh == nullptr)
    {
        OGSError::box("The VoxelGrid is faulty");  // write name of layer.
        return;
    }

    _mesh_model.addMesh(grid.release());
    this->done(QDialog::Accepted);
}

std::vector<std::string> Vtu2GridDialog::getSelectedObjects(QStringList list)
{
    std::vector<std::string> indexList;
    std::transform(list.begin(), list.end(), std::back_inserter(indexList),
                   [](auto const& index) { return index.toStdString(); });
    return indexList;
}