/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the TreeItem class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    : itemData_(std::move(data)), parentItem_(parent)
{
}

/**
 * A pointer to each of the child items belonging to this item
 * will be stored in the childItems_ private member variable. When
 * the class's destructor is called, it must delete each of these
 * to ensure that their memory is reused.
 */
TreeItem::~TreeItem()
{
    qDeleteAll(childItems_);
}

/**
 * Add a child to the tree item
 */
void TreeItem::appendChild(TreeItem* item)
{
    childItems_.append(item);
}

/**
 * Returns the child that corresponds to the specified row number
 * in the item's list of child items
 * Returns nullptr if that child does not exist.
 */
TreeItem* TreeItem::child(int row) const
{
    if (childItems_.count() > row)
    {
        return childItems_.value(row);
    }

    return nullptr;
}

/**
 * Returns the number of childItems_
 */
int TreeItem::childCount() const
{
    return childItems_.count();
}

/**
 * Returns the item's location within its parent's list of items.
 */
int TreeItem::row() const
{
    if (parentItem_)
    {
        return parentItem_->childItems_.indexOf(const_cast<TreeItem*>(this));
    }

    return 0;
}

/**
 * Returns the number of columns of data in the item.
 */
int TreeItem::columnCount() const
{
    return itemData_.count();
}

/**
 * Returns the data from all the columns.
 */
QVariant TreeItem::data(int column) const
{
    return itemData_.value(column);
}

/**
 * Sets the data at a given column.
 */
bool TreeItem::setData( int column, const QVariant &value )
{
    if (column < 0 || column >= itemData_.size())
    {
        return false;
    }

    itemData_[column] = value;
    return true;
}
/**
 * Returns the parent object of the tree item.
 */
TreeItem* TreeItem::parentItem() const
{
    return parentItem_;
}

/**
 * Removes a number of children of the TreeItem starting from the given position.
 */
bool TreeItem::removeChildren(int position, int count)
{
    if (position < 0 || position + count > childItems_.size())
    {
        return false;
    }

    for (int row = 0; row < count; ++row)
    {
        delete childItems_.takeAt(position);
    }

    return true;
}
