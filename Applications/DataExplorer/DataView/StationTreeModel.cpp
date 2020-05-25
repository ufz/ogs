/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the StationTreeModel class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "StationTreeModel.h"


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
    delete rootItem_;
    rootData << "Station Name" << "x" << "y" << "z";
    rootItem_ = new ModelTreeItem(rootData, nullptr, nullptr);
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
    {
        return QModelIndex();
    }

    ModelTreeItem* parentItem;

    if (!parent.isValid())
    {
        parentItem = static_cast<ModelTreeItem*>(rootItem_);
    }
    else
    {
        parentItem = static_cast<ModelTreeItem*>(parent.internalPointer());
    }

    auto* childItem = static_cast<ModelTreeItem*>(parentItem->child(row));
    if (childItem)
    {
        QModelIndex newIndex = createIndex(row, column, childItem);
        // assign ModelIndex to BaseItem so it can communicate with the model
        BaseItem* item = childItem->getItem();
        if (item != nullptr)
        {
            item->setModelIndex(newIndex);
        }
        return newIndex;
    }

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
        auto* treeItem = static_cast<ModelTreeItem*>(index.internalPointer());
        TreeItem* parentItem = treeItem->parentItem();
        listName = parentItem->data(0).toString();
        return treeItem->getStation();
    }

    return nullptr;
}

vtkPolyDataAlgorithm* StationTreeModel::vtkSource(const std::string &name) const
{
    std::size_t nLists = lists_.size();
    for (std::size_t i = 0; i < nLists; i++)
    {
        if (name == lists_[i]->data(0).toString().toStdString())
        {
            return dynamic_cast<BaseItem*>(lists_[i]->getItem())->vtkSource();
        }
    }
    return nullptr;
}

void StationTreeModel::setNameForItem(const std::string& stn_vec_name,
                                      std::size_t const id,
                                      std::string const& item_name)
{
    auto const stn_list = find_if(
        lists_.begin(), lists_.end(), [&stn_vec_name](ModelTreeItem* item) {
            return (stn_vec_name == item->data(0).toString().toStdString());
        });

    if (stn_list == lists_.end() ||
        id >= static_cast<std::size_t>((*stn_list)->childCount()))
    {
        return;
    }
    TreeItem *const item = (*stn_list)->child(id);
    item->setData(0, QString::fromStdString(item_name));
}

/**
 * Inserts a subtree under rootItem_.
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
    auto* group =
        new ModelTreeItem(grpName, rootItem_, new BaseItem(listName, stations));
    lists_.push_back(group);
    rootItem_->appendChild(group);
    int vectorSize = stations->size();

    for (int i = 0; i < vectorSize; i++)
    {
        QList<QVariant> stn;
        stn << QString::fromStdString(static_cast<GeoLib::Station*>((*stations)[i])->getName())
            << QString::number((*(*stations)[i])[0],'f')
            << QString::number((*(*stations)[i])[1],'f')
            << QString::number((*(*stations)[i])[2],'f');

        auto* child = new ModelTreeItem(stn, group);
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
        auto* item = static_cast<ModelTreeItem*>(getItem(index));

        // also delete the lists entry in the list directory of the model
        for (std::size_t i = 0; i < lists_.size(); i++)
        {
            if (item == lists_[i])
            {
                lists_.erase(lists_.begin() + i);
            }
        }

        removeRows(0, item->childCount(), index);
        removeRows(item->row(), 1, parent(index));
    }
}

/**
 * Removes the TreeItem with the given name including all its children
 */
void StationTreeModel::removeStationList(const std::string &name)
{
    for (auto& list : lists_)
    {
        if (name == list->data(0).toString().toStdString())
        {
            removeStationList(createIndex(list->row(), 0, list));
        }
    }
}
