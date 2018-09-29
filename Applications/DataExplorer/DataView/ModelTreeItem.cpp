/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Implementation of the ModelTreeItem class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ModelTreeItem.h"

ModelTreeItem::ModelTreeItem(const QList<QVariant> &data, TreeItem* parent, BaseItem* item)
    : TreeItem(data, parent), _item(item), _stn(nullptr)
{
}

BaseItem* ModelTreeItem::getItem() const
{
    if (_item != nullptr)
        return _item;
    return nullptr;
}

