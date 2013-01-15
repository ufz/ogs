/**
 * \file
 * \author Karsten Rink
 * \date   2011-11-22
 * \brief  Definition of the ProcessItem class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
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
