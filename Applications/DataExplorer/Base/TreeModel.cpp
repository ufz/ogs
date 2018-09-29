/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the TreeModel class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TreeModel.h"

#include "TreeItem.h"

#include <QModelIndex>
#include <QStringList>
#include <QVariant>

/**
 * A model for the QTreeView implementing the tree as a double-linked list.
 */
TreeModel::TreeModel( QObject* parent )
    : QAbstractItemModel(parent)
{
    //_modelType = TREE_MODEL;
    QList<QVariant> rootData;
    rootData << "1" << "2" << "3";
    _rootItem = new TreeItem(rootData, nullptr);
    //setupModelData(data, _rootItem);
}

TreeModel::~TreeModel()
{
    delete _rootItem;
}

/**
 * Returns the index of a TreeItem given its parent and position.
 * It is first checked if the model index is valid. If it is not, it is assumed that a
 * top-level item is being referred to; otherwise, the data pointer from the model
 * index is obtained with its internalPointer() function and is used to reference a
 * TreeItem object
 * \param row Row of the tree object
 * \param column Column of the tree object
 * \param parent Index of the parent object
 */
QModelIndex TreeModel::index(int row, int column, const QModelIndex &parent) const
{
    if (!hasIndex(row, column, parent))
        return QModelIndex();

    TreeItem* parentItem;

    if (!parent.isValid())
        parentItem = _rootItem;
    else
        parentItem = static_cast<TreeItem*>(parent.internalPointer());

    TreeItem* childItem = parentItem->child(row);
    if (childItem)
        return createIndex(row, column, childItem);

    return QModelIndex();
}

/**
 * Returns the model index of a TreeItem based on its index.
 */
QModelIndex TreeModel::parent(const QModelIndex &index) const
{
    if (!index.isValid())
        return QModelIndex();

    auto* childItem = static_cast<TreeItem*>(index.internalPointer());
    TreeItem* parentItem = childItem->parentItem();

    if (parentItem == _rootItem)
        return QModelIndex();

    return createIndex(parentItem->row(), 0, parentItem);
}

/**
 * Returns the number of child items for the TreeItem that corresponds to a given
 * model index, or the number of top-level items if an invalid index is specified.
 */
int TreeModel::rowCount(const QModelIndex &parent) const
{
    TreeItem* parentItem;
    if (parent.column() > 0)
        return 0;

    if (!parent.isValid())
        parentItem = _rootItem;
    else
        parentItem = static_cast<TreeItem*>(parent.internalPointer());

    return parentItem->childCount();
}

/**
 * Determines how many columns are present for a given model index.
 */
int TreeModel::columnCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return static_cast<TreeItem*>(parent.internalPointer())->columnCount();

    return _rootItem->columnCount();
}

void TreeModel::updateData()
{
}
/**
 * Since each item manages its own columns, the column number is used to retrieve
 * the data with the TreeItem::data() function
 */
QVariant TreeModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    if (role == Qt::EditRole || role == Qt::DisplayRole)
    {
        auto* item = static_cast<TreeItem*>(index.internalPointer());

        return item->data(index.column());
    }

    return QVariant();
}

bool TreeModel::setData( const QModelIndex &index, const QVariant &value, int role /* = Qt::EditRole */ )
{
    if (!index.isValid())
        return false;

    if (role == Qt::EditRole)
    {
        auto* item = static_cast<TreeItem*>(index.internalPointer());
        item->setData(index.column(), value);
        return true;
    }
    return false;
}
Qt::ItemFlags TreeModel::flags(const QModelIndex &index) const
{
    if (!index.isValid())
        return nullptr;

    return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}

/**
 * Returns the Item characterized by the given index.
 */
TreeItem* TreeModel::getItem(const QModelIndex &index) const
{
    if (index.isValid())
    {
        auto* item = static_cast<TreeItem*>(index.internalPointer());
        if (item)
            return item;
    }
    return _rootItem;
}

QVariant TreeModel::headerData(int section, Qt::Orientation orientation, int role) const
{
    if (orientation == Qt::Horizontal && role == Qt::DisplayRole)
        return _rootItem->data(section);

    return QVariant();
}

/**
 * Removes item from the model.
 */
bool TreeModel::removeRows(int position, int count, const QModelIndex & parent)
{
    TreeItem* parentItem = getItem(parent);
    bool success = true;

    beginRemoveRows(parent, position, position + count - 1);
    success = parentItem->removeChildren(position, count);
    endRemoveRows();

    return success;
}

/**
 * Test method for creating a tree model
 */
void TreeModel::setupModelData(const QStringList &lines, TreeItem* parent)
{
    QList<TreeItem*> parents;
    QList<int> indentations;
    parents << parent;
    indentations << 0;

    int number = 0;

    while (number < lines.count())
    {
        int position = 0;
        while (position < lines[number].length())
        {
            if (lines[number].mid(position, 1) != " ")
                break;
            position++;
        }

        QString lineData = lines[number].mid(position).trimmed();

        if (!lineData.isEmpty())
        {
            // Read the column data from the rest of the line.
            QStringList columnStrings = lineData.split("\t", QString::SkipEmptyParts);
            QList<QVariant> columnData;
            for (int column = 0; column < columnStrings.count(); ++column)
                columnData << columnStrings[column];

            if (position > indentations.last())
            {
                // The last child of the current parent is now the new parent
                // unless the current parent has no children.

                if (parents.last()->childCount() > 0)
                {
                    parents << parents.last()->child(parents.last()->childCount(
                                                             ) - 1);
                    indentations << position;
                }
            }
            else
                while (position < indentations.last() && parents.count() > 0)
                {
                    parents.pop_back();
                    indentations.pop_back();
                }

            // Append a new item to the current parent's list of children.
            parents.last()->appendChild(new TreeItem(columnData, parents.last()));
        }

        number++;
    }
}

TreeItem* TreeModel::rootItem() const
{
    return _rootItem;
}
