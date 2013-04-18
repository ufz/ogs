/**
 * \file
 * \author Karsten Rink
 * \date   2011-05-10
 * \brief  Implementation of the ElementTreeModel class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementTreeModel.h"
#include "OGSError.h"
#include "TreeItem.h"
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"

#include "VtkMeshSource.h"

/**
 * Constructor.
 */
ElementTreeModel::ElementTreeModel( QObject* parent )
	: TreeModel(parent), _mesh_source(nullptr)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Name" << "Type" << "" << "";
	_rootItem = new TreeItem(rootData, nullptr);
}

ElementTreeModel::~ElementTreeModel()
{
}

void ElementTreeModel::setElement(vtkUnstructuredGridAlgorithm const*const grid, const unsigned elem_index)
{
	this->_mesh_source = grid;
	this->clearView();

	VtkMeshSource const*const source = dynamic_cast<VtkMeshSource const*const>(grid);

	if (source)
	{
		const MeshLib::Mesh* mesh = source->GetMesh();
		const MeshLib::Element* elem = mesh->getElement(elem_index);

		QList<QVariant> elemData;
		elemData << "Element " + QString::number(elem_index) << "" << "" << "";
		TreeItem* elemItem = new TreeItem(elemData, _rootItem);
		_rootItem->appendChild(elemItem);

		QList<QVariant> typeData;
		typeData << "Element Type: " << QString::fromStdString(MshElemType2String(elem->getGeomType()));
		TreeItem* typeItem = new TreeItem(typeData, elemItem);
		elemItem->appendChild(typeItem);

		QList<QVariant> materialData;
		materialData << "MaterialID: " << QString::number(elem->getValue());
		TreeItem* matItem = new TreeItem(materialData, elemItem);
		elemItem->appendChild(matItem);

		QList<QVariant> volData;
		volData << "Area/Volume: " <<
		QString::number(mesh->getElement(elem_index)->getContent());
		TreeItem* volItem = new TreeItem(volData, elemItem);
		elemItem->appendChild(volItem);

		QList<QVariant> nodeListData;
		nodeListData << "Nodes" << "" << "" << "";
		TreeItem* nodeListItem = new TreeItem(nodeListData, elemItem);
		elemItem->appendChild(nodeListItem);

		//const std::vector<MeshLib::Node*> nodes_vec = grid->getNodes();
		size_t nElemNodes = elem->getNNodes();
		for (size_t i = 0; i < nElemNodes; i++)
		{
			const MeshLib::Node* node = elem->getNode(i);
			QList<QVariant> nodeData;
			nodeData << "Node " + QString::number(node->getID()) <<
			QString::number((*node)[0]) << QString::number((*node)[1]) <<
			QString::number((*node)[2]);
			TreeItem* nodeItem = new TreeItem(nodeData, nodeListItem);
			nodeListItem->appendChild(nodeItem);
		}
		reset();
	}
}

void ElementTreeModel::clearView()
{
	_rootItem->removeChildren(0, _rootItem->childCount());
	reset();
}

