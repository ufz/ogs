/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the ModelTreeItem class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ModelTreeItem.h"

ModelTreeItem::ModelTreeItem(const QList<QVariant> &data, TreeItem* parent, BaseItem* item)
    : TreeItem(data, parent), item_(item), stn_(nullptr)
{
}

BaseItem* ModelTreeItem::getItem() const
{
    if (item_ != nullptr)
    {
        return item_;
    }
    return nullptr;
}

