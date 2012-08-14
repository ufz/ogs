/**
 * \file ModelTreeItem.h
 * KR Initial implementation
 */

#ifndef QMODELTREEITEM_H
#define QMODELTREEITEM_H

#include "BaseItem.h"
#include "Station.h"
#include "TreeItem.h"

/**
 * \brief A TreeItem containing some additional information used in the StationModel.
 *
 * \sa TreeItem
 */
class ModelTreeItem : public TreeItem
{
public:
	/**
	 * Constructor
	 * \param data The data associated with each column
	 * \param parent The parent item in the tree
	 * \param item The ModelItem-object
	 */
	ModelTreeItem(const QList<QVariant> &data, TreeItem* parent, BaseItem* item = NULL);
	~ModelTreeItem() { delete _item; }

	/// Returns the station object from which this item has been constructed
	GEOLIB::Station* getStation() { return _stn; }

	/// Returns the BaseItem associated with this item
	BaseItem* getItem() const;

	/// Associates a station object with this item
	void setStation(GEOLIB::Station* stn) { _stn = stn; }

	/// Associates a BaseItem with this item
	void setItem( BaseItem* item ) { _item = item; }

private:
	BaseItem* _item;
	GEOLIB::Station* _stn;
};

#endif //QMODELTREEITEM_H
