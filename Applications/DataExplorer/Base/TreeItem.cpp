/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the TreeItem class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <utility>

#include "TreeItem.h"

/**
 * The constructor is only used to record the item's parent
 * and the data associated with each column.
 */
TreeItem::TreeItem(QList<QVariant> data, TreeItem* parent)
    : _itemData(std::move(data)), _parentItem(parent)
{
}

/**
 * A pointer to each of the child items belonging to this item
 * will be stored in the _childItems private member variable. When
 * the class's destructor is called, it must delete each of these
 * to ensure that their memory is reused.
 */
TreeItem::~TreeItem()
{
    qDeleteAll(_childItems);
}

/**
 * Add a child to the tree item
 */
void TreeItem::appendChild(TreeItem* item)
{
    _childItems.append(item);
}

/**
 * Returns the child that corresponds to the specified row number
 * in the item's list of child items
 * Returns nullptr if that child does not exist.
 */
TreeItem* TreeItem::child(int row) const
{
    if (_childItems.count() > row)
        return _childItems.value(row);
    else
        return nullptr;
}

/**
 * Returns the number of _childItems
 */
int TreeItem::childCount() const
{
    return _childItems.count();
}

/**
 * Returns the item's location within its parent's list of items.
 */
int TreeItem::row() const
{
    if (_parentItem)
        return _parentItem->_childItems.indexOf(const_cast<TreeItem*>(this));

    return 0;
}

/**
 * Returns the number of columns of data in the item.
 */
int TreeItem::columnCount() const
{
    return _itemData.count();
}

/**
 * Returns the data from all the columns.
 */
QVariant TreeItem::data(int column) const
{
    return _itemData.value(column);
}

/**
 * Sets the data at a given column.
 */
bool TreeItem::setData( int column, const QVariant &value )
{
    if (column < 0 || column >= _itemData.size())
        return false;

    _itemData[column] = value;
    return true;
}
/**
 * Returns the parent object of the tree item.
 */
TreeItem* TreeItem::parentItem() const
{
    return _parentItem;
}

/**
 * Removes a number of children of the TreeItem starting from the given position.
 */
bool TreeItem::removeChildren(int position, int count)
{
    if (position < 0 || position + count > _childItems.size())
        return false;

    for (int row = 0; row < count; ++row)
        delete _childItems.takeAt(position);

    return true;
}
