/**
 * \file CondItem.h
 * 20/10/2010 KR Initial implementation
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
	const GEOLIB::GeoObject* getGeoObject() { return this->getGeoObject(); }

private:
	const FEMCondition* _item;
};

#endif //CONDITEM_H
