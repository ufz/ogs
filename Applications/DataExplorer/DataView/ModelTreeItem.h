/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the ModelTreeItem class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "GeoLib/Station.h"

#include "BaseItem.h"
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
    ModelTreeItem(const QList<QVariant>& data, TreeItem* parent,
                  BaseItem* item = nullptr);
    ~ModelTreeItem() override { delete item_; }
    /// Returns the station object from which this item has been constructed
    GeoLib::Station* getStation() { return stn_; }

    /// Returns the BaseItem associated with this item
    BaseItem* getItem() const;

    /// Associates a station object with this item
    void setStation(GeoLib::Station* stn) { stn_ = stn; }

    /// Associates a BaseItem with this item
    void setItem( BaseItem* item ) { item_ = item; }

private:
    BaseItem* item_;
    GeoLib::Station* stn_;
};
