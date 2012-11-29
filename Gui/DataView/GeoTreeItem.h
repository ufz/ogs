/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file GeoTreeItem.h
 *
 * Created on 2010-12-08 by Karsten Rink
 */

#ifndef GEOTREEITEM_H
#define GEOTREEITEM_H

#include "TreeItem.h"

#include "GeoObject.h"

/**
 * \brief A TreeItem containing an additional GeoObject
 *
 * \sa TreeItem
 */
class GeoTreeItem : public TreeItem
{
public:
	/**
	 * Constructor
	 * \param data The data associated with each column
	 * \param parent The parent item in the tree
	 * \param item The GeoObject (i.e. Point, Polyline or Surface)
	 */
	GeoTreeItem(const QList<QVariant> &data,
	            TreeItem* parent,
	            const GeoLib::GeoObject* item = NULL) : TreeItem(data, parent), _item(item) {}
	~GeoTreeItem() {}

	/// Returns the geo-object associated with this item (i.e. a point, polyline or surface).
	const GeoLib::GeoObject* getGeoObject() const { return _item; }

private:
	const GeoLib::GeoObject* _item;
};

#endif //GEOTREEITEM_H
