/**
 * \file
 * \date   2023-05-11
 * \brief  Implementation of the AddFaultsToVoxelGridDialog class.
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AddFaultsToVoxelGridDialog.h"

#include <QStringList>
#include <QStringListModel>
#include <algorithm>
#include <string>

#include "Base/OGSError.h"
#include "Base/Utils.h"
#include "MeshLib/Mesh.h"
#include "MeshModel.h"
#include "MeshToolsLib/MeshGenerators/AddFaultToVoxelGrid.h"

QString no_voxel_str = "[No voxel grid available.]";

AddFaultsToVoxelGridDialog::AddFaultsToVoxelGridDialog(MeshModel& mesh_model,
                                                       QDialog* parent)
    : QDialog(parent), _mesh_model(mesh_model)
{
    setupUi(this);
    QStringList voxelGridList;
    QStringList faultsList;

    for (int model_index = 0; model_index < mesh_model.rowCount();
         ++model_index)
    {
        auto const* mesh = mesh_model.getMesh(mesh_model.index(model_index, 0));
        assert(mesh);
        if (MeshToolsLib::MeshGenerator::AddFaultToVoxelGrid::isVoxelGrid(
                *mesh))
        {
            voxelGridList.append(QString::fromStdString(mesh->getName()));
        }
        if (mesh->getDimension() == 2)
        {
            faultsList.append(QString::fromStdString(mesh->getName()));
        }
    }
    if (voxelGridList.empty())
    {
        voxelGridList.append(no_voxel_str);
    }

    if (faultsList.empty())
    {
        faultsList.append("[No 2D faults available.]");
    }

    _voxelGrids.setStringList(voxelGridList);
    _meshes2D.setStringList(faultsList);
    this->voxelGridListBox->addItems(_voxelGrids.stringList());
    this->all2Dmeshes->setModel(&_meshes2D);
    this->selectedFaults->setModel(&_selFaults);
}

void AddFaultsToVoxelGridDialog::on_selectDataButton_pressed()
{
    Utils::moveSelectedItems(this->all2Dmeshes, _meshes2D, _selFaults);
}

void AddFaultsToVoxelGridDialog::on_deselectDataButton_pressed()
{
    Utils::moveSelectedItems(this->selectedFaults, _selFaults, _meshes2D);
}

void AddFaultsToVoxelGridDialog::accept()
{
    using namespace MeshToolsLib::MeshGenerator;
    if (this->voxelGridListBox->currentText() == no_voxel_str)
    {
        OGSError::box("Please specify the input voxel grid.");
        return;
    }

    auto* input_voxelgrid = _mesh_model.getMesh(
        this->voxelGridListBox->currentText().toStdString());
    assert(input_voxelgrid);
    std::unique_ptr<MeshLib::Mesh> voxelgrid(
        new MeshLib::Mesh(*input_voxelgrid));
    assert(voxelgrid);

    auto const* mat_ids = MeshLib::materialIDs(*voxelgrid);

    if (mat_ids == nullptr)
    {
        OGSError::box("Input mesh has no material IDs");
        return;
    }

    std::vector<std::string> const selected_faults =
        Utils::getSelectedObjects(_selFaults.stringList());

    if (this->_selFaults.rowCount() == 0)
    {
        OGSError::box("Please specify the fault(s) to be added to the mesh.");
        return;
    }

    auto max_mat_id = std::max_element(mat_ids->cbegin(), mat_ids->cend());
    assert(max_mat_id);
    auto fault_id = *max_mat_id;
    for (auto const& fault_name : selected_faults)
    {
        fault_id++;
        auto const* fault = _mesh_model.getMesh(fault_name);
        if (AddFaultToVoxelGrid::addFaultToVoxelGrid(voxelgrid.get(), fault,
                                                     fault_id))
        {
            INFO("The fault '{}' was added.", fault_name);
            continue;
        }
        OGSError::box("The fault " + QString::fromStdString(fault_name) +
                      " could not be added.");
    }

    _mesh_model.addMesh(voxelgrid.release());
    this->done(QDialog::Accepted);
}