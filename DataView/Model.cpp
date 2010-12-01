/**
 * \file Model.cpp
 * 24/9/2009 LB Initial implementation
 * 05/05/2010 KR 2d graphic functionality removed and various layout changes
 *
 * Implementation of Model
 */

// ** INCLUDES **
#include "Model.h"
#include "TreeItem.h"

#include <vtkPolyDataAlgorithm.h>

#include <QItemSelection>
#include <QDebug>

#include <iostream>


Model::Model( QString name, QObject* parent /*= 0*/ )
: TreeModel(parent), _name(name), _vtkSource(NULL)
{
}
/*
ModelItem* Model::itemFromIndex( const QModelIndex& index ) const
{
	if (index.isValid())
		return (ModelItem*)(index.internalPointer());
	else
		return NULL;
}

QModelIndex Model::indexFromItem( const ModelItem* item ) const
{
	if (item != NULL)
	{
		for (int i = 0; i < rowCount(); ++i)
		{
			QModelIndex modelIndex = index(i, 0);
			ModelItem* itemFromModel = itemFromIndex(modelIndex);

			if (item == itemFromModel)
				return modelIndex;
		}
	}
	return QModelIndex();
}
*/
QVector<Model*> Model::subModels() const
{
	return _subModels;
}

Qt::ItemFlags Model::flags( const QModelIndex& index ) const
{
	if (!index.isValid())
		return Qt::ItemIsEnabled;

	return QAbstractItemModel::flags(index) | Qt::ItemIsEditable;
}

bool Model::removeRows( int row, int count, const QModelIndex & parent /*= QModelIndex() */ )
{
	QItemSelection deletedRowsSelection;
	beginInsertColumns(parent, row, row + count - 1);
	for (int i = 0; i < count; i++)
	{
		//QVector<ModelItem*>::iterator it = _data.begin() + row;
		//ModelItem* item = *it;
		//QModelIndex index = indexFromItem(item);
		//QModelIndex lastIndex = index.sibling(index.row(), columnCount()-1);
		//deletedRowsSelection.select(index, lastIndex);
		//delete item;
		//_data.erase(it);
	}
 	_selectedItems.merge(deletedRowsSelection, QItemSelectionModel::Deselect);
	endInsertColumns();

	return true;
}

void Model::updateData()
{
	foreach( Model* subModel, _subModels )
		subModel->updateData();
	reset();

	//_vtkSource->Update();
}

void Model::clearData()
{
	_rootItem->removeChildren(0, _rootItem->childCount());
	emit selectionCleared();
}

void Model::clearSelectedData()
{
	while (_selectedItems.indexes().size() > 0)
	{
		QModelIndexList indices = _selectedItems.indexes();
		int lastRow = 0;
		for (int i = 0; i < indices.size(); i++)
			if (indices.at(i).row() > lastRow)
				lastRow = indices.at(i).row();

		qDebug() << "Model: Removing row" << lastRow;
		removeRow(lastRow);
	}

	//_vtkSource->Update();

	emit selectionCleared();
}
