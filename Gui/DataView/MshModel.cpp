/**
 * \file
 * \author Lars Bilke
 * \date   2009-10-19
 * \brief  Implementation of the MshModel class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MshModel.h"

// ThirdParty/logog
#include "logog/include/logog.hpp"

// ** INCLUDES **
#include "MshItem.h"
#include "StringTools.h"
#include "TreeItem.h"
#include "VtkMeshSource.h"

// MeshLib
#include "Node.h"
#include "Elements/Element.h"
#include "MeshEnums.h"

// Qt
#include <QFileInfo>
#include <QString>

#include "StringTools.h"

MshModel::MshModel(ProjectData &project, QObject* parent /*= 0*/ )
	: TreeModel(parent), _project(project)
{
	delete _rootItem;
	QList<QVariant> rootData;
	rootData << "Mesh Name" << "Type" << "Node IDs";
	_rootItem = new TreeItem(rootData, nullptr);
}

int MshModel::columnCount( const QModelIndex &parent /*= QModelIndex()*/ ) const
{
	Q_UNUSED(parent)

	return 3;
}

void MshModel::addMesh(MeshLib::Mesh* mesh)
{
	_project.addMesh(mesh);
	this->addMeshObject(mesh);
}

void MshModel::addMeshObject(const MeshLib::Mesh* mesh)
{
	INFO("name: %s", mesh->getName().c_str());
	QString display_name (QString::fromStdString(mesh->getName()));
	QList<QVariant> meshData;
	meshData << display_name << "";
	MshItem* newMesh = new MshItem(meshData, _rootItem, mesh);
	if (newMesh->vtkSource())
		newMesh->vtkSource()->SetName(display_name);
	_rootItem->appendChild(newMesh);

	// display elements
	const std::vector<MeshLib::Element*> &elems = mesh->getElements();
	const size_t nElems (elems.size());
	QString elem_type_string("");
	MeshElemType elem_type(MeshElemType::INVALID);

	for (size_t i = 0; i < nElems; i++)
	{
		const MeshLib::Element* current_element (elems[i]);
		MeshElemType t (current_element->getGeomType());
		QList<QVariant> elemData;
		elemData.reserve(3);
		if (t != elem_type)
		{
			elem_type = t;
			elem_type_string = QString::fromStdString(MeshElemType2String(t));
		}

		QString nodestxt(QString::number(current_element->getNode(0)->getID()));
		const size_t nNodes(current_element->getNNodes());
		for (size_t j = 1; j < nNodes; j++)
			nodestxt.append(", " + QString::number(current_element->getNode(j)->getID()));

		elemData << "Element " + QString::number(i) << elem_type_string << nodestxt;
		newMesh->appendChild(new TreeItem(elemData, newMesh));
	}

	reset();
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
			emit meshRemoved(this, this->index(i, 0, QModelIndex()));
			_rootItem->removeChildren(i,1);
			reset();
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
	const std::vector<MeshLib::Mesh*> msh_vec = _project.getMeshObjects();
	for (std::vector<MeshLib::Mesh*>::const_iterator it(msh_vec.begin()); it != msh_vec.end(); ++it)
		if (!this->getMesh((*it)->getName())) // if Mesh is not yet added to GUI, do it now
			addMeshObject(*it);
}

VtkMeshSource* MshModel::vtkSource(const QModelIndex &idx) const
{
	if (idx.isValid())
	{
		MshItem* item = static_cast<MshItem*>(this->getItem(idx));
		return item->vtkSource();
	}

	INFO("MshModel::vtkSource(): Specified index does not exist.");
	return nullptr;
}

VtkMeshSource* MshModel::vtkSource(const std::string &name) const
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

