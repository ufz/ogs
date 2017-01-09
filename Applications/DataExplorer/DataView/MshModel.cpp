/**
 * \file
 * \author Lars Bilke
 * \date   2009-10-19
 * \brief  Implementation of the MshModel class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MshModel.h"

// Qt
#include <QFileInfo>
#include <QString>

#include <vtkUnstructuredGridAlgorithm.h>

#include <logog/include/logog.hpp>

#include "MshItem.h"
#include "StringTools.h"
#include "TreeItem.h"

// MeshLib
#include "MeshLib/Node.h"
#include "Elements/Element.h"

const QVariant MshModel::element_str = "Element";
const std::map<MeshLib::MeshElemType, QVariant> MshModel::elem_type_map = MshModel::createMeshElemTypeMap();

MshModel::MshModel(DataHolderLib::Project &project, QObject* parent /*= 0*/ )
    : TreeModel(parent), _project(project)
{
    delete _rootItem;
    QList<QVariant> rootData;
    rootData << "Mesh Name" << "#" << "Type";
    _rootItem = new TreeItem(rootData, nullptr);
}

int MshModel::columnCount( const QModelIndex &parent /*= QModelIndex()*/ ) const
{
    Q_UNUSED(parent)

    return 3;
}

void MshModel::addMesh(std::unique_ptr<MeshLib::Mesh> mesh)
{
    _project.addMesh(std::move(mesh));
    auto const& meshes(_project.getMeshObjects());
    this->addMeshObject(meshes.back().get());
}

void MshModel::addMesh(MeshLib::Mesh* mesh)
{
    _project.addMesh(std::unique_ptr<MeshLib::Mesh>(mesh));
    auto const& meshes(_project.getMeshObjects());
    this->addMeshObject(meshes.back().get());
}

void MshModel::addMeshObject(const MeshLib::Mesh* mesh)
{
    beginResetModel();

    INFO("name: %s", mesh->getName().c_str());
    QVariant const display_name (QString::fromStdString(mesh->getName()));
    QList<QVariant> meshData;
    meshData << display_name << "" << "";
    MshItem *const newMesh = new MshItem(meshData, _rootItem, mesh);
    _rootItem->appendChild(newMesh);

    // display elements
    std::vector<MeshLib::Element*> const& elems = mesh->getElements();
    int const nElems (static_cast<int>(elems.size()));

    for (int i = 0; i < nElems; i++)
    {
        QList<QVariant> elemData;
        elemData.reserve(3);
        elemData << element_str << i << elem_type_map.at(elems[i]->getGeomType());
        newMesh->appendChild(new TreeItem(elemData, newMesh));
    }

    endResetModel();
    emit meshAdded(this, this->index(_rootItem->childCount() - 1, 0, QModelIndex()));
}

const MeshLib::Mesh* MshModel::getMesh(const QModelIndex &idx) const
{
    if (idx.isValid())
    {
        MshItem* item = dynamic_cast<MshItem*>(this->getItem(idx));
        if (item)
            return item->getMesh();
        else
            return nullptr;
    }
    WARN("MshModel::getMesh(): Specified index does not exist.");
    return nullptr;
}

const MeshLib::Mesh* MshModel::getMesh(const std::string &name) const
{
    for (int i = 0; i < _rootItem->childCount(); i++)
    {
        MshItem* item = static_cast<MshItem*>(_rootItem->child(i));
        if (item->data(0).toString().toStdString().compare(name) == 0)
            return item->getMesh();
    }

    INFO("MshModel::getMesh(): No entry found with name \"%s\".", name.c_str());
    return nullptr;
}

bool MshModel::removeMesh(const QModelIndex &idx)
{
    if (idx.isValid())
    {
        MshItem* item = dynamic_cast<MshItem*>(this->getItem(idx));
        if (item)
            return this->removeMesh(item->getMesh()->getName());
        return false;
    }
    return false;
}

bool MshModel::removeMesh(const std::string &name)
{
    for (int i = 0; i < _rootItem->childCount(); i++)
    {
        TreeItem* item = _rootItem->child(i);
        if (item->data(0).toString().toStdString().compare(name) == 0)
        {
            beginResetModel();
            emit meshRemoved(this, this->index(i, 0, QModelIndex()));
            _rootItem->removeChildren(i,1);
            endResetModel();
            return _project.removeMesh(name);
        }
    }

    INFO("MshModel::removeMesh(): No entry found with name \"%s\".", name.c_str());
    return false;
}

void MshModel::updateMesh(MeshLib::Mesh* mesh)
{
    for (int i = 0; i < _rootItem->childCount(); i++)
    {
        if (dynamic_cast<MshItem*>(this->_rootItem->child(i))->getMesh() == mesh)
        {
            emit meshRemoved(this, this->index(i, 0, QModelIndex()));
            _rootItem->removeChildren(i,1);
        }
    }
    this->addMeshObject(mesh);
}

void MshModel::updateModel()
{
    auto const& mesh_vec = _project.getMeshObjects();
    for (auto const& mesh : mesh_vec)
        if (!getMesh(mesh->getName())) // if Mesh is not yet added to GUI, do it now
            addMeshObject(mesh.get());
}

std::map<MeshLib::MeshElemType, QVariant> MshModel::createMeshElemTypeMap()
{
    std::vector<MeshLib::MeshElemType> const& elem_types (MeshLib::getMeshElemTypes());
    std::map<MeshLib::MeshElemType, QVariant> elem_map;

    for (MeshLib::MeshElemType t : elem_types)
        elem_map[t] = QVariant(QString::fromStdString(MeshLib::MeshElemType2String(t)));

    return elem_map;
}

vtkUnstructuredGridAlgorithm* MshModel::vtkSource(const QModelIndex &idx) const
{
    if (idx.isValid())
    {
        MshItem* item = static_cast<MshItem*>(this->getItem(idx));
        return item->vtkSource();
    }

    INFO("MshModel::vtkSource(): Specified index does not exist.");
    return nullptr;
}

vtkUnstructuredGridAlgorithm* MshModel::vtkSource(const std::string &name) const
{
    for (int i = 0; i < _rootItem->childCount(); i++)
    {
        MshItem* item = static_cast<MshItem*>(_rootItem->child(i));
        if (item->data(0).toString().toStdString().compare(name) == 0)
            return item->vtkSource();
    }

    INFO("MshModel::vtkSource(): No entry found with name \"%s\".", name.c_str());
    return nullptr;
}
