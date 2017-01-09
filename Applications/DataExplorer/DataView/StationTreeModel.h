/**
 * \file
 * \author Karsten Rink
 * \date   no date
 * \brief  Definition of the StationTreeModel class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "ModelTreeItem.h"
#include "Point.h"
#include "TreeModel.h"
#include <vtkPolyDataAlgorithm.h>

namespace GeoLib
{
class Station;
class StationBorehole;
}

class QString;
class QModelIndex;
class PropertyBounds;

/**
 * \brief A model for the StationTreeView implementing a tree as a double-linked list.
 *
 * A model for the StationTreeView implementing a tree as a double-linked list.
 * In addition to a simple TreeModel each item also contains a 2D / 3D GraphicsItem for visualization.
 * \sa TreeModel, StationTreeView, TreeItem, ModelTreeItem
 */
class StationTreeModel : public TreeModel
{
    Q_OBJECT

public:
    StationTreeModel(QObject* parent = 0);
    ~StationTreeModel();

    void addStationList(QString listName, const std::vector<GeoLib::Point*>* stations);
    const std::vector<ModelTreeItem*> &getLists() { return _lists; }
    QModelIndex index(int row, int column, const QModelIndex &parent = QModelIndex()) const;
    //BaseItem* itemFromIndex( const QModelIndex& index ) const;
    void removeStationList(QModelIndex index);
    void removeStationList(const std::string &name);
    void setNameForItem(const std::string& stn_vec_name, std::size_t id,
                        std::string const& item_name);
    GeoLib::Station* stationFromIndex( const QModelIndex& index, QString &listName ) const;
    vtkPolyDataAlgorithm* vtkSource(const std::string &name) const;

private:
    std::vector<ModelTreeItem*> _lists;
};
