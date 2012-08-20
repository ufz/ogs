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

// MeshLib
#include "Node.h"
#include "Elements/Element.h"
#include "MshEnums.h"

// Qt
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

void MshModel::addMesh(MeshLib::Mesh* mesh)
{
	_project.addMesh(mesh);
	this->addMeshObject(mesh);
}

void MshModel::addMeshObject(MeshLib::Mesh* mesh)
{
	std::cout << "name: " << mesh->getName() << std::endl;
	QString display_name (QString::fromStdString(mesh->getName()));
	QList<QVariant> meshData;
	meshData << display_name << "";
	MshItem* newMesh = new MshItem(meshData, _rootItem, mesh);
	if (newMesh->vtkSource())
		newMesh->vtkSource()->SetName(display_name);
	_rootItem->appendChild(newMesh);

	// display elements
	const std::vector<MeshLib::Element*> elems = mesh->getElements();
	const size_t nElems (elems.size());
	QString elem_type_string("");
	MshElemType::type elem_type(MshElemType::INVALID);

	for (size_t i = 0; i < nElems; i++)
	{
		const MeshLib::Element* current_element (elems[i]);
		MshElemType::type t (current_element->getType());
		QList<QVariant> elemData;
		if (t != elem_type)
		{
			elem_type = t;
			elem_type_string = QString::fromStdString(MshElemType2String(t));
		}

		QString nodestxt("");
		const size_t nNodes(current_element->getNNodes());
		for (size_t j = 0; j < nNodes; j++)
			nodestxt.append(QString::number(current_element->getNode(j)->getID()) + ", ");

		elemData << "Element " + QString::number(i) << elem_type_string << nodestxt.left(nodestxt.length() - 2);

		TreeItem* elem = new TreeItem(elemData, newMesh);
		newMesh->appendChild(elem);
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
			return NULL;
	}
	std::cout << "MshModel::getMesh() - Specified index does not exist." << std::endl;
	return NULL;
}

const MeshLib::Mesh* MshModel::getMesh(const std::string &name) const
{
	for (int i = 0; i < _rootItem->childCount(); i++)
	{
		MshItem* item = static_cast<MshItem*>(_rootItem->child(i));
		if (item->data(0).toString().toStdString().compare(name) == 0)
			return item->getMesh();
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


