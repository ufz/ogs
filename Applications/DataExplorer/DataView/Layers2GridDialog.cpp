/**
 * \file
 * \date   2023-04-26
 * \brief  Implementation of the Layers2GridDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Layers2GridDialog.h"

#include <QStringList>
#include <QStringListModel>

#include "Base/StrictDoubleValidator.h"
#include "Base/Utils.h"
#include "GeoLib/AABB.h"
#include "MeshLib/IO/writeMeshToFile.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshModel.h"
#include "MeshToolsLib/MeshGenerators/VoxelGridFromLayeredMeshes.h"

Layers2GridDialog::Layers2GridDialog(MeshModel& mesh_model, QDialog* parent)
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

    _layeredMeshes.setStringList(MeshList);
    this->allMeshView->setModel(&_layeredMeshes);
    this->allMeshView->setDragDropMode(QAbstractItemView::InternalMove);
}

void Layers2GridDialog::on_deleteMeshButton_pressed()
{
    QModelIndex const selected =
        this->allMeshView->selectionModel()->currentIndex();
    _layeredMeshes.removeRow(selected.row());
    QStringList list = _layeredMeshes.stringList();
    _layeredMeshes.setStringList(list);
}

void Layers2GridDialog::on_orderButton_pressed()
{
    _layeredMeshes.sort(0);
}

void Layers2GridDialog::on_upOrderButton_pressed()
{
    QModelIndex selected = this->allMeshView->selectionModel()->currentIndex();
    QStringList list = _layeredMeshes.stringList();
    int row = selected.row();
    if (row > 0)
    {
        QString list_item = list[row - 1];
        list[row - 1] = selected.data().toString();
        list[row] = list_item;
    }
    _layeredMeshes.setStringList(list);
    this->allMeshView->selectionModel()->setCurrentIndex(
        _layeredMeshes.index(row - 1), QItemSelectionModel::SelectCurrent);
}

void Layers2GridDialog::on_downOrderButton_pressed()
{
    QModelIndex selected = this->allMeshView->selectionModel()->currentIndex();
    QStringList list = _layeredMeshes.stringList();
    int row = selected.row();
    if (row < list.size() - 1 && row != -1)
    {
        QString list_item = list[row + 1];
        list[row + 1] = selected.data().toString();
        list[row] = list_item;
    }
    _layeredMeshes.setStringList(list);
    this->allMeshView->selectionModel()->setCurrentIndex(
        _layeredMeshes.index(row + 1), QItemSelectionModel::SelectCurrent);
}

void Layers2GridDialog::updateExpectedVoxel()
{
    QString const xin = this->xlineEdit->text();
    QString const yin = this->ylineEdit->text();
    QString const zin = this->zlineEdit->text();
    bool ok;
    double const xinput = xin.toDouble();
    double const yinput = (yin.toDouble(&ok)) ? yin.toDouble() : xin.toDouble();
    double const zinput = (zin.toDouble(&ok)) ? zin.toDouble() : xin.toDouble();

    if (_layeredMeshes.stringList()[0] == "[No Mesh available.]")
    {
        this->expectedVoxelLabel->setText("approximated Voxel: undefined");
        return;
    }
    if (xin.isEmpty() || xinput == 0)
    {
        this->expectedVoxelLabel->setText("approximated Voxel: undefined");
        return;
    }

    std::vector<std::string> layered_meshes =
        Utils::getSelectedObjects(_layeredMeshes.stringList());
    auto* const mesh_top = _mesh_model.getMesh(layered_meshes.front());
    auto* const mesh_bottom = _mesh_model.getMesh(layered_meshes.back());
    auto const& nodes_top = mesh_top->getNodes();

    GeoLib::AABB const aabb_top(nodes_top.cbegin(), nodes_top.cend());
    auto const& nodes_bottom = mesh_bottom->getNodes();
    GeoLib::AABB const aabb_bottom(nodes_bottom.cbegin(), nodes_bottom.cend());
    auto const min_b = aabb_bottom.getMinPoint();
    auto const max_b = aabb_bottom.getMaxPoint();
    auto const max_t = aabb_top.getMaxPoint();
    double const expectedVoxel = (max_b[0] - min_b[0]) * (max_b[1] - min_b[1]) *
                                 (max_t[2] - min_b[2]) / xinput / yinput /
                                 zinput;

    int const exponent = std::floor(std::log10(abs(expectedVoxel)));
    this->expectedVoxelLabel->setText(
        "approximated Voxel = " +
        QString::number(std::round(expectedVoxel / std::pow(10, exponent))) +
        " x 10^" + QString::number(exponent));
}

void Layers2GridDialog::on_xlineEdit_textChanged()
{
    updateExpectedVoxel();
}

void Layers2GridDialog::on_ylineEdit_textChanged()
{
    updateExpectedVoxel();
}

void Layers2GridDialog::on_zlineEdit_textChanged()
{
    updateExpectedVoxel();
}

void Layers2GridDialog::accept()
{
    if (this->_layeredMeshes.rowCount() == 1)
    {
        OGSError::box(
            "Please specify the input meshes. \n At least two layers are "
            "required to create a 3D Mesh");
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

    std::vector<std::string> layered_meshes =
        Utils::getSelectedObjects(_layeredMeshes.stringList());

    std::vector<const MeshLib::Mesh*> layers;
    layers.reserve(layered_meshes.size());

    for (auto const& layer : layered_meshes)
    {
        auto mesh(_mesh_model.getMesh(layer));
        if (mesh == nullptr)
        {
            OGSError::box("Input layer " + QString::fromStdString(layer) +
                          " not found. Aborting...");
            return;
        }
        layers.push_back(mesh);
    }

    bool const dilate = this->dilateBox->isChecked();
    std::array<double, 3> const cellsize = {xinput, yinput, zinput};
    constexpr double minval = std::numeric_limits<double>::max();
    constexpr double maxval = std::numeric_limits<double>::lowest();
    std::pair<MathLib::Point3d, MathLib::Point3d> extent(
        MathLib::Point3d{{minval, minval, minval}},
        MathLib::Point3d{{maxval, maxval, maxval}});
    auto mesh(MeshToolsLib::MeshGenerators::VoxelFromLayeredMeshes::
                  createVoxelFromLayeredMesh(extent, layers, cellsize, dilate));

    if (mesh == nullptr)
    {
        OGSError::box("The VoxelGrid is faulty");
        return;
    }
    OGSError::box("The VoxelGrid is fine");

    _mesh_model.addMesh(mesh.release());
    this->done(QDialog::Accepted);
}
