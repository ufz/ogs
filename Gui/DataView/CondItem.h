/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file CondItem.h
 *
 * Created on 2010-10-20 by Karsten Rink
 */

#ifndef CONDITEM_H
#define CONDITEM_H

#include "FEMCondition.h"
#include "TreeItem.h"
#include "VtkPointsSource.h"

/**
 * \brief A TreeItem containing a condition of a FEM (BC, IC or ST).
 * \sa TreeItem
 */
class CondItem : public TreeItem
{
public:
	/// Constructor
	CondItem(const QList<QVariant> &data, TreeItem* parent, const FEMCondition* cond)
		: TreeItem(data, parent), _item(cond)
	{
	}

	~CondItem() {}

	/// Returns the FEM Condition associated with the item.
	const FEMCondition* getItem() { return _item; }

	/// Returns the geo-object on which the condition is placed.
	const GeoLib::GeoObject* getGeoObject() { return this->getGeoObject(); }

private:
	const FEMCondition* _item;
};

#endif //CONDITEM_H
