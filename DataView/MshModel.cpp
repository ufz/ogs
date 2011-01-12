/**
 * \file MshModel.cpp
 * 19/10/2009 LB Initial implementation
 * 12/05/2010 KR re-implementation
 *
 * Implementation of MshModel
 */

// ** INCLUDES **
#include "MshModel.h"
#include "MshItem.h"
#include "StringTools.h"
#include "TreeItem.h"
#include "VtkMeshSource.h"

#include "msh_lib.h"

#include <QString>
#include <QFileInfo>
#include <QTime>

MshModel::MshModel(ProjectData &project, QObject* parent /*= 0*/ )
: TreeModel(parent), _project(project)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Mesh Name";
	_rootItem = new TreeItem(rootData, NULL);
}

void MshModel::addMesh(Mesh_Group::CFEMesh* mesh, std::string &name)
{
	_project.addMesh(mesh, name);
	this->addMesh(new GridAdapter(mesh), name);
}

void MshModel::addMesh(GridAdapter* mesh, std::string &name)
{
	std::cout << "name: " << name << std::endl;
	QFileInfo fi(QString::fromStdString(name));
	name = fi.baseName().toStdString();
	mesh->setName(name);
	QList<QVariant> meshData;
	meshData.push_back(QVariant(QString::fromStdString(name)));
	MshItem* newMesh = new MshItem(meshData, _rootItem, mesh);
	if (newMesh->vtkSource())
		newMesh->vtkSource()->SetName(fi.fileName());
	std::cout << "name: " << fi.fileName().toStdString() << std::endl;
	_rootItem->appendChild(newMesh);
	reset();

	emit meshAdded(this, this->index(_rootItem->childCount()-1, 0, QModelIndex()));
}

const GridAdapter* MshModel::getMesh(const QModelIndex &idx) const
{
	if (idx.isValid())
	{
		MshItem* item = static_cast<MshItem*>(this->getItem(idx));
		return item->getGrid();
	}
	std::cout << "MshModel::getMesh() - Specified index does not exist." << std::endl;
	return NULL;
}

const GridAdapter* MshModel::getMesh(const std::string &name) const
{
	for (int i=0; i<_rootItem->childCount(); i++)
	{
		MshItem* item = static_cast<MshItem*>(_rootItem->child(i));
		if (item->data(0).toString().toStdString().compare(name) == 0)
			return item->getGrid();
	}

	std::cout << "MshModel::getMesh() - No entry found with name \"" << name << "\"." << std::endl;
	return NULL;
}


bool MshModel::removeMesh(const QModelIndex &idx)
{
	if (idx.isValid())
	{
		MshItem* item = static_cast<MshItem*>(this->getItem(idx));
		emit meshRemoved(this, idx);
		_rootItem->removeChildren(item->row(),1);
		reset();
		return true;
	}

	std::cout << "MshModel::removeMesh() - Specified index does not exist." << std::endl;
	return false;
}

bool MshModel::removeMesh(const std::string &name)
{
	for (int i=0; i<_rootItem->childCount(); i++)
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

	std::cout << "MshModel::removeMesh() - No entry found with name \"" << name << "." << std::endl;
	return false;
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
	for (int i=0; i<_rootItem->childCount(); i++)
	{
		MshItem* item = static_cast<MshItem*>(_rootItem->child(i));
		if (item->data(0).toString().toStdString().compare(name) == 0)
			return item->vtkSource();
	}

	std::cout << "MshModel::getMesh() - No entry found with name \"" << name << "\"." << std::endl;
	return NULL;
}

int MshModel::columnCount( const QModelIndex& parent /*= QModelIndex()*/ ) const
{
	Q_UNUSED(parent)

	return 1;
}
/*
bool MshModel::isUniqueMeshName(std::string &name)
{
	int count(0);
	bool isUnique(false);
	std::string cpName;

	while (!isUnique)
	{
		isUnique = true;
		cpName = name;

		count++;
		// If the original name already exists we start to add numbers to name for
		// as long as it takes to make the name unique.
		if (count>1) cpName = cpName + "-" + number2str(count);

		for (int i=0; i<_rootItem->childCount(); i++)
		{
			TreeItem* item = _rootItem->child(i);
			if (item->data(0).toString().toStdString().compare(cpName) == 0) isUnique = false;
		}
	}

	// At this point cpName is a unique name and isUnique is true.
	// If cpName is not the original name, "name" is changed and isUnique is set to false,
	// indicating that a vector with the original name already exists.
	if (count>1)
	{
		isUnique = false;
		name = cpName;
	}
	return isUnique;
}
*/

Mesh_Group::CFEMesh* MshModel::loadMeshFromFile(std::string fileName)
{
	std::cout << "FEMRead ... " << std::flush;
#ifndef NDEBUG
	QTime myTimer;
	myTimer.start();
#endif
	FEMDeleteAll();

	CFEMesh* msh = FEMRead(fileName.substr(0, fileName.length()-4));
	if (msh)
	{
#ifndef NDEBUG
		QTime constructTimer;
		constructTimer.start();
#endif
		msh->ConstructGrid();
		std::cout << "Nr. Nodes: " << msh->nod_vector.size() << endl;
		std::cout << "Nr. Elements: " << msh->ele_vector.size() << endl;

#ifndef NDEBUG
		std::cout << "constructGrid time: " << constructTimer.elapsed() << " ms" << std::endl;
#endif
		msh->FillTransformMatrix();

#ifndef NDEBUG
		std::cout << "Loading time: " << myTimer.elapsed() << " ms" << std::endl;
#endif
		return msh;
	}

    cout << "Failed to load a mesh file: " << fileName << endl;
	return NULL;
}
