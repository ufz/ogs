/**
 * \file ModelTreeItem.cpp
 * KR Initial implementation
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

