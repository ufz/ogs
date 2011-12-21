/**
 * \file MshModel.cpp
 * 19/10/2009 LB Initial implementation
 * 12/05/2010 KR re-implementation
 *
 * Implementation of MshModel
 */

// ** INCLUDES **
#include "MshItem.h"
#include "MshModel.h"
#include "StringTools.h"
#include "TreeItem.h"
#include "VtkMeshSource.h"

#include "msh_lib.h"

#include <QFileInfo>
#include <QString>
#include <QTime>

MshModel::MshModel(ProjectData &project, QObject* parent /*= 0*/ )
	: TreeModel(parent), _project(project)
{
	delete _rootItem;
	QList<QVariant> rootData;
	rootData << "Mesh Name" << "Type" << "Node IDs";
	_rootItem = new TreeItem(rootData, NULL);
}

int MshModel::columnCount( const QModelIndex &parent /*= QModelIndex()*/ ) const
{
	Q_UNUSED(parent)

	return 3;
}

void MshModel::addMesh(MeshLib::CFEMesh* mesh, std::string &name)
{
	_project.addMesh(mesh, name);
	GridAdapter* grid = new GridAdapter(mesh);
	grid->setName(name);
	this->addMeshObject(grid);
}

void MshModel::addMesh(GridAdapter* mesh)
{
	MeshLib::CFEMesh* ogsmesh( const_cast<MeshLib::CFEMesh*>(mesh->getCFEMesh()) );
	std::string msh_name = mesh->getName();
	_project.addMesh(ogsmesh, msh_name);
	this->addMeshObject(mesh);
}

void MshModel::addMeshObject(GridAdapter* mesh)
{
	std::string name(mesh->getName());
	std::cout << "name: " << name << std::endl;
	QFileInfo fi(QString::fromStdString(name));
	name = fi.baseName().toStdString();
	mesh->setName(name);
	QList<QVariant> meshData;
	meshData << QString::fromStdString(name) << "";
	MshItem* newMesh = new MshItem(meshData, _rootItem, mesh);
	if (newMesh->vtkSource())
		newMesh->vtkSource()->SetName(fi.fileName());
	_rootItem->appendChild(newMesh);

	// display elements
	const std::vector<GridAdapter::Element*>* elems = mesh->getElements();
	size_t nElems (elems->size());
	for (size_t i = 0; i < nElems; i++)
	{
		QList<QVariant> elemData;
		elemData << "Element " + QString::number(i) << QString::fromStdString(
		        MshElemType2String((*elems)[i]->type));

		QString nodestxt("");
		size_t nNodes((*elems)[i]->nodes.size());
		for (size_t j = 0; j < nNodes; j++)
			nodestxt.append(QString::number((*elems)[i]->nodes[j]) + ", ");
		elemData << nodestxt.left(nodestxt.length() - 2);

		TreeItem* elem = new TreeItem(elemData, newMesh);
		newMesh->appendChild(elem);
	}

	reset();

	emit meshAdded(this, this->index(_rootItem->childCount() - 1, 0, QModelIndex()));
}

const GridAdapter* MshModel::getMesh(const QModelIndex &idx) const
{
	if (idx.isValid())
	{
		MshItem* item = dynamic_cast<MshItem*>(this->getItem(idx));
		if (item)
			return item->getGrid();
		else
			return NULL;
	}
	std::cout << "MshModel::getMesh() - Specified index does not exist." << std::endl;
	return NULL;
}

const GridAdapter* MshModel::getMesh(const std::string &name) const
{
	for (int i = 0; i < _rootItem->childCount(); i++)
	{
		MshItem* item = static_cast<MshItem*>(_rootItem->child(i));
		if (item->data(0).toString().toStdString().compare(name) == 0)
			return item->getGrid();
	}

	std::cout << "MshModel::getMesh() - No entry found with name \"" << name << "\"." <<
	std::endl;
	return NULL;
}

bool MshModel::removeMesh(const QModelIndex &idx)
{
	if (idx.isValid())
	{
		MshItem* item = dynamic_cast<MshItem*>(this->getItem(idx));
		if (item)
		{
			emit meshRemoved(this, idx);
			_rootItem->removeChildren(item->row(),1);
			reset();
			return true;
		}
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

	std::cout << "MshModel::removeMesh() - No entry found with name \"" << name << "." <<
	std::endl;
	return false;
}

void MshModel::updateModel()
{
	const std::map<std::string, MeshLib::CFEMesh*> msh_vec = _project.getMeshObjects();
	for (std::map<std::string, MeshLib::CFEMesh*>::const_iterator it(msh_vec.begin());
	     it != msh_vec.end(); ++it)
		if (this->getMesh(it->first) == NULL) // if GridAdapter does not yet exist, create one.
		{
			std::string name = it->first;
			addMeshObject(new GridAdapter(it->second));
		}
}

VtkMeshSource* MshModel::vtkSource(const QModelIndex &idx) const
{
	if (idx.isValid())
	{
		MshItem* item = static_cast<MshItem*>(this->getItem(idx));
		return item->vtkSource();
	}

	std::cout << "MshModel::removeMesh() - Specified index does not exist." << std::endl;
	return NULL;
}

VtkMeshSource* MshModel::vtkSource(const std::string &name) const
{
	for (int i = 0; i < _rootItem->childCount(); i++)
	{
		MshItem* item = static_cast<MshItem*>(_rootItem->child(i));
		if (item->data(0).toString().toStdString().compare(name) == 0)
			return item->vtkSource();
	}

	std::cout << "MshModel::getMesh() - No entry found with name \"" << name << "\"." <<
	std::endl;
	return NULL;
}


