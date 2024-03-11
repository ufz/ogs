/**
 * \file
 * \date   2023-05-11
 * \brief  Implementation of the Vtu2GridDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Vtu2GridDialog.h"

#include <vtkXMLUnstructuredGridWriter.h>

#include <QStringList>
#include <QStringListModel>
#include <optional>
#include <string>

#include "Base/OGSError.h"
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
        this->expectedVoxelLabel->setText(
            "Expected number of voxels: undefined");
    }

    _allMeshes.setStringList(MeshList);
    this->meshListBox->addItems(_allMeshes.stringList());
    this->xlineEdit->setFocus();
}

std::optional<std::array<double, 3>> fillXYZ(QString xin, QString yin,
                                             QString zin)
{
    bool ok;
    if (!xin.toDouble(&ok))
    {
        return std::nullopt;
    }
    double const xinput = xin.toDouble();
    double const yinput = (yin.toDouble(&ok)) ? yin.toDouble() : xinput;
    double const zinput = (zin.toDouble(&ok)) ? zin.toDouble() : xinput;

    if (xinput <= 0 || yinput <= 0 || zinput <= 0)
    {
        return std::nullopt;
    }

    return std::optional<std::array<double, 3>>{{xinput, yinput, zinput}};
}

Eigen::Vector3d getMeshExtent(MeshLib::Mesh const* _mesh)
{
    auto const& nodes = _mesh->getNodes();
    GeoLib::AABB const aabb(nodes.cbegin(), nodes.cend());
    auto const& min = aabb.getMinPoint();
    auto const& max = aabb.getMaxPoint();
    return max - min;
}

void Vtu2GridDialog::updateExpectedVoxel()
{
    if (_allMeshes.stringList()[0] == "[No Mesh available.]")
    {
        this->expectedVoxelLabel->setText("approximated Voxel: undefined");
        return;
    }

    auto const opt_xyz =
        fillXYZ(this->xlineEdit->text(), this->ylineEdit->text(),
                this->zlineEdit->text());

    if (!opt_xyz)
    {
        this->expectedVoxelLabel->setText("approximated Voxel: undefined");
        return;
    }

    auto const delta = getMeshExtent(
        _mesh_model.getMesh(this->meshListBox->currentText().toStdString()));
    double const expected_voxel = (delta[0]) * (delta[1]) * (delta[2]) /
                                  (*opt_xyz)[0] / (*opt_xyz)[1] / (*opt_xyz)[2];

    int const exponent = std::floor(std::log10(std::abs(expected_voxel)));
    this->expectedVoxelLabel->setText(
        "approximated Voxel = " +
        QString::number(std::round(expected_voxel / std::pow(10, exponent))) +
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

    auto cellsize = fillXYZ(this->xlineEdit->text(), this->ylineEdit->text(),
                            this->zlineEdit->text());

    if (!cellsize)
    {
        OGSError::box(
            "At least the x-length of a voxel must be specified and > 0.\n If "
            "y-/z-input "
            "are not specified, equal to 0, or not a real number, they are "
            "treated as "
            "the x-input.");
    }
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
        VoxelGridFromMesh::getNumberOfVoxelPerDimension(ranges, *cellsize);
    std::unique_ptr<MeshLib::Mesh> grid(
        generateRegularHexMesh(dims[0], dims[1], dims[2], (*cellsize)[0],
                               (*cellsize)[1], (*cellsize)[2], min, "grid"));
    if (grid == nullptr)
    {
        OGSError::box(QString::fromStdString(
            fmt::format("Could not generate regular hex mesh. With "
                        "parameters dims={} {} {}, cellsize={} {} {}",
                        dims[0], dims[1], dims[2], (*cellsize)[0],
                        (*cellsize)[1], (*cellsize)[2])));
    }

    std::vector<int>* cell_ids =
        grid->getProperties().createNewPropertyVector<int>(
            VoxelGridFromMesh::cell_id_name, MeshLib::MeshItemType::Cell, 1);
    if (cell_ids == nullptr)
    {
        OGSError::box("Could not create cell ids.");
    }
    *cell_ids = VoxelGridFromMesh::assignCellIds(mesh, min, dims, *cellsize);
    if (!VoxelGridFromMesh::removeUnusedGridCells(mesh, grid))
    {
        return;
    }

    VoxelGridFromMesh::mapMeshArraysOntoGrid(mesh, grid);

    if (grid == nullptr)
    {
        OGSError::box("No voxelgrid could be created from the mesh.");
        return;
    }

    _mesh_model.addMesh(grid.release());
    this->done(QDialog::Accepted);
}
