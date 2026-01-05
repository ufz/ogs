// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ModelTreeItem.h"

ModelTreeItem::ModelTreeItem(const QList<QVariant>& data, TreeItem* parent,
                             BaseItem* item)
    : TreeItem(data, parent), _item(item), _stn(nullptr)
{
}

BaseItem* ModelTreeItem::getItem() const
{
    if (_item != nullptr)
    {
        return _item;
    }
    return nullptr;
}
