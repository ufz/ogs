/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file ModelTreeItem.cpp
 *
 * Created on by Karsten Rink
 */

#include "ModelTreeItem.h"

ModelTreeItem::ModelTreeItem(const QList<QVariant> &data, TreeItem* parent, BaseItem* item)
	: TreeItem(data, parent), _item(item)
{
}

BaseItem* ModelTreeItem::getItem() const
{
	if (_item != NULL)
		return _item;
	return NULL;
}

