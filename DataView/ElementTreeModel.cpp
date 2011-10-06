/**
 * \file ElementTreeModel.cpp
 * 2011/05/10 KR Initial implementation
 */

#include "ElementTreeModel.h"
#include "OGSError.h"
#include "TreeItem.h"

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

void ElementTreeModel::setElement(const GridAdapter* grid, const size_t elem_index)
{
	this->clearView();
	GridAdapter::Element* elem = (*(grid->getElements()))[elem_index];

	QList<QVariant> elemData;
	elemData << "Element " + QString::number(elem_index) << "" << "" << "";
	TreeItem* elemItem = new TreeItem(elemData, _rootItem);
	_rootItem->appendChild(elemItem);

	QList<QVariant> typeData;
	typeData << "Element Type: " << QString::fromStdString(MshElemType2String(elem->type));
	TreeItem* typeItem = new TreeItem(typeData, elemItem);
	elemItem->appendChild(typeItem);

	QList<QVariant> materialData;
	materialData << "MaterialID: " << QString::number(elem->material);
	TreeItem* matItem = new TreeItem(materialData, elemItem);
	elemItem->appendChild(matItem);

	QList<QVariant> volData;
	volData << "Area/Volume: " <<
	QString::number(grid->getCFEMesh()->getElementVector()[elem_index]->calcVolume());
	TreeItem* volItem = new TreeItem(volData, elemItem);
	elemItem->appendChild(volItem);

	QList<QVariant> nodeListData;
	nodeListData << "Nodes" << "" << "" << "";
	TreeItem* nodeListItem = new TreeItem(nodeListData, elemItem);
	elemItem->appendChild(nodeListItem);

	const std::vector<GEOLIB::Point*>* nodes = grid->getNodes();
	for (size_t i = 0; i < elem->nodes.size(); i++)
	{
		QList<QVariant> nodeData;
		nodeData << "Node " + QString::number(elem->nodes[i]) <<
		QString::number((*(*nodes)[i])[0]) << QString::number((*(*nodes)[i])[1]) <<
		QString::number((*(*nodes)[i])[2]);
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

