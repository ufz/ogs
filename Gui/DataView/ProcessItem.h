/**
 * \file ProcessItem.h
 * 2011/11/22 KR Initial implementation
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
