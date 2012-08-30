/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file StationTreeModel.cpp
 *
 * Created on by Karsten Rink
 */

#include "BaseItem.h"
#include "OGSError.h"
#include "Station.h"
#include "StationTreeModel.h"

#include <QDebug>

/**
 * Constructor.
 */
StationTreeModel::StationTreeModel( QObject* parent )
	: TreeModel(parent)
{
	QList<QVariant> rootData;
	delete _rootItem;
	rootData << "Station Name" << "x" << "y" << "z";
	_rootItem = new ModelTreeItem(rootData, NULL, NULL);
}

StationTreeModel::~StationTreeModel()
{
}

/**
 * Returns the model index of an item in the tree.
 * \param row The row where the item is located
 * \param column The column where the item is located
 * \param parent The parent of the item
 * \return The model index of the item
 */
QModelIndex StationTreeModel::index( int row, int column,
                                     const QModelIndex &parent /*= QModelIndex()*/ ) const
{
	if (!hasIndex(row, column, parent))
		return QModelIndex();

	ModelTreeItem* parentItem;

	if (!parent.isValid())
		parentItem = (ModelTreeItem*)(_rootItem);
	else
		parentItem = static_cast<ModelTreeItem*>(parent.internalPointer());

	ModelTreeItem* childItem = (ModelTreeItem*)(parentItem->child(row));
	if (childItem)
	{
		QModelIndex newIndex = createIndex(row, column, childItem);
		// assign ModelIndex to BaseItem so it can communicate with the model
		BaseItem* item = childItem->getItem();
		if ( item != NULL )
			item->setModelIndex(newIndex);
		return newIndex;
	}
	else
		return QModelIndex();
}

/*
 * Returns the model item associated with a tree item. This is a re-implementation of
 * QOgsModelBase::itemFromIndex to deal with TreeItems with an optionally attached BaseItem
 * \param index Index of the requested item
 * \return The BaseItem associated with the tree item
 *
   BaseItem* StationTreeModel::itemFromIndex( const QModelIndex& index ) const
   {
    if (index.isValid())
    {
        ModelTreeItem* treeItem = static_cast<ModelTreeItem*>(index.internalPointer());
        BaseItem* baseItem = treeItem->getItem();
        return baseItem;
    }
    else
        return NULL;
   }
 */

/**
 * Returns the Station-Object of the ModelTreeItem with the given index and the name of the list this station belongs to.
 * \param index Index of the requested item
 * \param listName Here, the method will put the name of the list this station belongs to.
 * \return The station object associated with the tree item
 */
GeoLib::Station* StationTreeModel::stationFromIndex( const QModelIndex& index,
                                                     QString &listName ) const
{
	if (index.isValid())
	{
		ModelTreeItem* treeItem = static_cast<ModelTreeItem*>(index.internalPointer());
		TreeItem* parentItem = treeItem->parentItem();
		listName = parentItem->data(0).toString();
		return treeItem->getStation();
	}
	else
		return NULL;
}

vtkPolyDataAlgorithm* StationTreeModel::vtkSource(const std::string &name) const
{
	size_t nLists = _lists.size();
	for (size_t i = 0; i < nLists; i++)
		if ( name.compare( _lists[i]->data(0).toString().toStdString() ) == 0 )
			return dynamic_cast<BaseItem*>(_lists[i]->getItem())->vtkSource();
	return NULL;
}

/**
 * Inserts a subtree under _rootItem.
 * \param listName Name of the new subtree. If no name is given a default name is assigned.
 * \param stations The list with stations to be added as children of that subtree
 */
void StationTreeModel::addStationList(QString listName, const std::vector<GeoLib::Point*>* stations)
{
	QList<QVariant> grpName;
	if (listName.compare("") == 0) // if no name is given a default name is assigned
	{
		listName = "List";
		listName.append(QString::number(rowCount() + 1));
	}
	grpName.push_back(QVariant(listName));
	grpName.push_back(QVariant(""));
	grpName.push_back(QVariant(""));
	grpName.push_back(QVariant(""));
	BaseItem* grpItem = new BaseItem(listName, stations);
	ModelTreeItem* group = new ModelTreeItem(grpName, _rootItem, grpItem);
	_lists.push_back(group);
	_rootItem->appendChild(group);
	int vectorSize = stations->size();

	for (int i = 0; i < vectorSize; i++)
	{
		QList<QVariant> stn;
		stn.push_back(QVariant(QString::fromStdString(static_cast<GeoLib::Station*>((*
		                                                                             stations)
		                                                                            [i])->
		                                              getName())));
		stn.push_back(QVariant(QString::number((*(*stations)[i])[0],'f')));
		stn.push_back(QVariant(QString::number((*(*stations)[i])[1],'f')));
		stn.push_back(QVariant(QString::number((*(*stations)[i])[2],'f')));

		ModelTreeItem* child = new ModelTreeItem(stn, group);
		child->setStation(static_cast<GeoLib::Station*>((*stations)[i]));
		group->appendChild(child);
	}

	qDebug() << "List" << listName << "loaded, " << stations->size() << "items added.";

	reset();
}

/**
 * Removes the TreeItem with the given Index including all its children
 */
void StationTreeModel::removeStationList(QModelIndex index)
{
	if (index.isValid()) //
	{
		ModelTreeItem* item = static_cast<ModelTreeItem*>(getItem(index));

		// also delete the lists entry in the list directory of the model
		for (size_t i = 0; i < _lists.size(); i++)
			if (item == _lists[i])
				_lists.erase(_lists.begin() + i);

		removeRows(0, item->childCount(), index);
		removeRows(item->row(), 1, parent(index));
	}
}

/**
 * Removes the TreeItem with the given name including all its children
 */
void StationTreeModel::removeStationList(const std::string &name)
{
	for (size_t i = 0; i < _lists.size(); i++)
		if ( name.compare( _lists[i]->data(0).toString().toStdString() ) == 0 )
			removeStationList(createIndex(_lists[i]->row(), 0, _lists[i]));
}

/*
 * Traverses the (sub-)tree of the TreeItem with the given index and looks for an item in that subtree
 * with the given name. Using this method is probably rather slow and not recommended.
 * If this method is called with an empty QModelIndex, the TreeModel will assign the rootItem to that
 * index.
 * \param index The TreeItem whose subtree will be searched.
 * \param name The name of the item that should be found
 * \return The QModelIndex of the desired item (or an empty index if no such item exists).
 *
   QModelIndex StationTreeModel::getItemByName(const QModelIndex &idx, const std::string &name) const
   {
    QModelIndex didx;
    TreeItem* item = getItem(idx);

    int nChildren = item->childCount();
    for (int i=0; i<nChildren; i++)
    {
        TreeItem* child = item->child(i);
        QString test = child->data(0).toString();

        if ( name.compare(child->data(0).toString().toStdString()) != 0 )
            didx = getItemByName( index(i, 0, idx), name );
        else didx = index(i, 0, idx);
        if (didx.isValid())
            return didx;

    }

    return QModelIndex(); // this is no valid QModelIndex and signifies that no item by the given name could be found.
   }
 */

/**
 * Filters the station list based on the property boundaries given in bounds.
 * Technically, the complete station list is removed from the model and only those items are re-loaded that fit the description.
 * If no station in the list fulfills the given description an error msg is given.
 */
void StationTreeModel::filterStations(const std::string &listName,
                                      const std::vector<GeoLib::Point*>* stations,
                                      const std::vector<PropertyBounds> &bounds)
{
	std::vector<GeoLib::Point*>* filteredStations = new std::vector<GeoLib::Point*>;

	size_t vectorSize = stations->size();
	for (size_t i = 0; i < vectorSize; i++)
		if (static_cast<GeoLib::Station*>((*stations)[i])->inSelection(bounds))
			filteredStations->push_back((*stations)[i]);

	if (filteredStations->empty())
		OGSError::box("No object is within the given boundaries.");  //The filtered list is empty.
	else
	{
		removeStationList(listName);
		this->addStationList(QString::fromStdString(listName), filteredStations);
		std::cout << "Filter applied to List \"" << listName << "\", " <<
		filteredStations->size() << " items added.";
	}
}
