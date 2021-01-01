/**
 * \file
 * \author Lars Bilke
 * \date   2009-10-19
 * \brief  Implementation of the MeshModel class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MeshModel.h"

#include <vtkUnstructuredGridAlgorithm.h>
#include <QFileInfo>
#include <QString>
#include "BaseLib/Logging.h"

#include "Applications/DataHolderLib/Project.h"
#include "BaseLib/StringTools.h"
#include "Elements/Element.h"
#include "MeshLib/Node.h"

#include "MeshItem.h"
#include "TreeItem.h"


const QVariant MeshModel::element_str = "Element";
const std::map<MeshLib::MeshElemType, QVariant> MeshModel::elem_type_map = MeshModel::createMeshElemTypeMap();

MeshModel::MeshModel(DataHolderLib::Project &project, QObject* parent /*= 0*/ )
    : TreeModel(parent), _project(project)
{
    delete _rootItem;
    QList<QVariant> rootData;
    rootData << "Mesh Name" << "#" << "Type";
    _rootItem = new TreeItem(rootData, nullptr);
}

int MeshModel::columnCount( const QModelIndex &parent /*= QModelIndex()*/ ) const
{
    Q_UNUSED(parent)

    return 3;
}

void MeshModel::addMesh(std::unique_ptr<MeshLib::Mesh> mesh)
{
    _project.addMesh(std::move(mesh));
    auto const& meshes(_project.getMeshObjects());
    this->addMeshObject(meshes.back().get());
}

void MeshModel::addMesh(MeshLib::Mesh* mesh)
{
    _project.addMesh(std::unique_ptr<MeshLib::Mesh>(mesh));
    auto const& meshes(_project.getMeshObjects());
    this->addMeshObject(meshes.back().get());
}

void MeshModel::addMeshObject(const MeshLib::Mesh* mesh)
{
    beginResetModel();

    INFO("name: {:s}", mesh->getName());
    QVariant const display_name (QString::fromStdString(mesh->getName()));
    QList<QVariant> meshData;
    meshData << display_name << "" << "";
    auto* const newMesh = new MeshItem(meshData, _rootItem, mesh);
    _rootItem->appendChild(newMesh);

    // display elements
    std::vector<MeshLib::Element*> const& elems = mesh->getElements();
    auto const nElems(static_cast<int>(elems.size()));

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

const MeshLib::Mesh* MeshModel::getMesh(const QModelIndex &idx) const
{
    if (idx.isValid())
    {
        auto* item = dynamic_cast<MeshItem*>(this->getItem(idx));
        if (item)
        {
            return item->getMesh();
        }

        return nullptr;
    }
    WARN("MeshModel::getMesh(): Specified index does not exist.");
    return nullptr;
}

const MeshLib::Mesh* MeshModel::getMesh(const std::string &name) const
{
    for (int i = 0; i < _rootItem->childCount(); i++)
    {
        auto* item = static_cast<MeshItem*>(_rootItem->child(i));
        if (item->data(0).toString().toStdString() == name)
        {
            return item->getMesh();
        }
    }

    INFO("MeshModel::getMesh(): No entry found with name \"{:s}\".", name);
    return nullptr;
}

bool MeshModel::removeMesh(const QModelIndex &idx)
{
    if (idx.isValid())
    {
        auto* item = dynamic_cast<MeshItem*>(this->getItem(idx));
        if (item)
        {
            return this->removeMesh(item->getMesh()->getName());
        }
        return false;
    }
    return false;
}

bool MeshModel::removeMesh(const std::string &name)
{
    for (int i = 0; i < _rootItem->childCount(); i++)
    {
        TreeItem* item = _rootItem->child(i);
        if (item->data(0).toString().toStdString() == name)
        {
            beginResetModel();
            emit meshRemoved(this, this->index(i, 0, QModelIndex()));
            _rootItem->removeChildren(i,1);
            endResetModel();
            return _project.removeMesh(name);
        }
    }

    INFO("MeshModel::removeMesh(): No entry found with name \"{:s}\".", name);
    return false;
}

void MeshModel::updateMesh(MeshLib::Mesh* mesh)
{
    for (int i = 0; i < _rootItem->childCount(); i++)
    {
        if (dynamic_cast<MeshItem*>(this->_rootItem->child(i))->getMesh() == mesh)
        {
            emit meshRemoved(this, this->index(i, 0, QModelIndex()));
            _rootItem->removeChildren(i,1);
        }
    }
    this->addMeshObject(mesh);
}

void MeshModel::updateModel()
{
    auto const& mesh_vec = _project.getMeshObjects();
    for (auto const& mesh : mesh_vec)
    {
        if (!getMesh(mesh->getName()))
        {  // if Mesh is not yet added to GUI, do it now
            addMeshObject(mesh.get());
        }
    }
}

std::map<MeshLib::MeshElemType, QVariant> MeshModel::createMeshElemTypeMap()
{
    std::vector<MeshLib::MeshElemType> const& elem_types (MeshLib::getMeshElemTypes());
    std::map<MeshLib::MeshElemType, QVariant> elem_map;

    for (MeshLib::MeshElemType t : elem_types)
    {
        elem_map[t] =
            QVariant(QString::fromStdString(MeshLib::MeshElemType2String(t)));
    }

    return elem_map;
}

vtkUnstructuredGridAlgorithm* MeshModel::vtkSource(const QModelIndex &idx) const
{
    if (idx.isValid())
    {
        auto* item = static_cast<MeshItem*>(this->getItem(idx));
        return item->vtkSource();
    }

    INFO("MeshModel::vtkSource(): Specified index does not exist.");
    return nullptr;
}

vtkUnstructuredGridAlgorithm* MeshModel::vtkSource(const std::string &name) const
{
    for (int i = 0; i < _rootItem->childCount(); i++)
    {
        auto* item = static_cast<MeshItem*>(_rootItem->child(i));
        if (item->data(0).toString().toStdString() == name)
        {
            return item->vtkSource();
        }
    }

    INFO("MeshModel::vtkSource(): No entry found with name \"{:s}\".", name);
    return nullptr;
}
