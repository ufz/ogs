/**
 * \file
 * \author Lars Bilke
 * \date   2010-06-23
 * \brief  Implementation of the TreeModelIterator class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license

 */

// ** INCLUDES **
#include "TreeModelIterator.h"

#include "TreeItem.h"
#include "TreeModel.h"

TreeModelIterator::TreeModelIterator(TreeModel* model)
    : current_(nullptr), model_(model)
{
    if (model_->rootItem()->childCount() > 0)
    {
        current_ = model_->rootItem();
    }
}

TreeItem* TreeModelIterator::operator*() const
{
    return current_;
}

TreeModelIterator& TreeModelIterator::operator++()
{
    if (current_)
    {
        current_ = next(current_);
    }

    return *this;
}

TreeItem* TreeModelIterator::next( const TreeItem* current )
{
    if (!current)
    {
        return nullptr;
    }

    TreeItem* next = nullptr;
    if (current->childCount())
    {
        // walk the child
        parentIndex_.push(currentIndex_);
        currentIndex_ = 0;
        next = current->child(0);
    }
    else
    {
        // walk the sibling
        TreeItem* parent = current->parentItem();
        next = parent ? parent->child(currentIndex_ + 1)
               : model_->rootItem()->child(currentIndex_ + 1);
        while (!next && parent)
        {
            // if we had no sibling walk up the parent
            // and try the sibling of that
            parent = parent->parentItem();
            currentIndex_ = parentIndex_.pop();
            next = parent ? parent->child(currentIndex_ + 1)
                   : model_->rootItem()->child(currentIndex_ + 1);
        }
        if (next)
        {
            ++(currentIndex_);
        }
    }
    return next;
}
