/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file ProcessItem.h
 *
 * Created on 2011-11-22 by Karsten Rink
 */

#ifndef PROCESSITEM_H
#define PROCESSITEM_H

#include "TreeItem.h"
#include "ProcessInfo.h"

/**
 * \brief A TreeItem representing process information.
 * \sa TreeItem
 */
class ProcessItem : public TreeItem
{
public:
	/// Constructor
	ProcessItem(const QList<QVariant> &data, TreeItem* parent, const ProcessInfo* pcs)
		: TreeItem(data, parent), _item(pcs)
	{
	}

	~ProcessItem() {}

	/// Returns the	Process Information associated with the item.
	const ProcessInfo* getItem() { return _item; }

private:
	const ProcessInfo* _item;
};

#endif //PROCESSITEM_H
