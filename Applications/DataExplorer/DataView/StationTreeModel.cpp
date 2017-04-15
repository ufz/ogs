/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the StationTreeModel class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "StationTreeModel.h"

#include <logog/include/logog.hpp>

#include "BaseItem.h"
#include "OGSError.h"
#include "Station.h"

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
    _rootItem = new ModelTreeItem(rootData, nullptr, nullptr);
}

StationTreeModel::~StationTreeModel() = default;

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
        if ( item != nullptr )
            item->setModelIndex(newIndex);
        return newIndex;
    }
    else
        return QModelIndex();
}

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
        return nullptr;
}

vtkPolyDataAlgorithm* StationTreeModel::vtkSource(const std::string &name) const
{
    std::size_t nLists = _lists.size();
    for (std::size_t i = 0; i < nLists; i++)
        if ( name.compare( _lists[i]->data(0).toString().toStdString() ) == 0 )
            return dynamic_cast<BaseItem*>(_lists[i]->getItem())->vtkSource();
    return nullptr;
}

void StationTreeModel::setNameForItem(const std::string& stn_vec_name,
                                      std::size_t id,
                                      std::string const& item_name)
{
    auto const stn_list = find_if(
        _lists.begin(), _lists.end(), [&stn_vec_name](ModelTreeItem* item) {
            return (stn_vec_name == item->data(0).toString().toStdString());
        });

    if (id >=  (*stn_list)->childCount())
        return;
    TreeItem *const item = (*stn_list)->child(id);
    item->setData(0, QString::fromStdString(item_name));
}

/**
 * Inserts a subtree under _rootItem.
 * \param listName Name of the new subtree. If no name is given a default name is assigned.
 * \param stations The list with stations to be added as children of that subtree
 */
void StationTreeModel::addStationList(QString listName, const std::vector<GeoLib::Point*>* stations)
{
    beginResetModel();

    QList<QVariant> grpName;
    if (listName.compare("") == 0) // if no name is given a default name is assigned
    {
        listName = "List";
        listName.append(QString::number(rowCount() + 1));
    }
    grpName << listName << "" << "" << "";
    ModelTreeItem* group = new ModelTreeItem(grpName, _rootItem, new BaseItem(listName, stations));
    _lists.push_back(group);
    _rootItem->appendChild(group);
    int vectorSize = stations->size();

    for (int i = 0; i < vectorSize; i++)
    {
        QList<QVariant> stn;
        stn << QString::fromStdString(static_cast<GeoLib::Station*>((*stations)[i])->getName())
            << QString::number((*(*stations)[i])[0],'f')
            << QString::number((*(*stations)[i])[1],'f')
            << QString::number((*(*stations)[i])[2],'f');

        ModelTreeItem* child = new ModelTreeItem(stn, group);
        child->setStation(static_cast<GeoLib::Station*>((*stations)[i]));
        group->appendChild(child);
    }

    qDebug() << "List" << listName << "loaded, " << stations->size() << "items added.";

    endResetModel();
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
        for (std::size_t i = 0; i < _lists.size(); i++)
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
    for (std::size_t i = 0; i < _lists.size(); i++)
        if ( name.compare( _lists[i]->data(0).toString().toStdString() ) == 0 )
            removeStationList(createIndex(_lists[i]->row(), 0, _lists[i]));
}
