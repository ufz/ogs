/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file ElementTreeModel.cpp
 *
 * Created on 2011-05-10 by Karsten Rink
 */

#include "ElementTreeModel.h"
#include "OGSError.h"
#include "TreeItem.h"
#include "Mesh.h"
#include "Node.h"
#include "Elements/Element.h"

/**
 * Constructor.
 */
ElementTreeModel::ElementTreeModel( QObject* parent )
	: TreeModel(parent)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Name" << "Type" << "" << "";
	_rootItem = new TreeItem(rootData, NULL);
}

ElementTreeModel::~ElementTreeModel()
{
}

void ElementTreeModel::setElement(const MeshLib::Mesh* grid, const size_t elem_index)
{
	this->clearView();
	const MeshLib::Element* elem = grid->getElement(elem_index);

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
	QString::number(grid->getElement(elem_index)->getContent());
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

void ElementTreeModel::clearView()
{
	_rootItem->removeChildren(0, _rootItem->childCount());
	reset();
}

