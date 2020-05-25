/**
 * \file
 * \author Karsten Rink
 * \date   2010-01-20
 * \brief  Definition of the BaseItem class.
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "VtkStationSource.h"
#include <QModelIndex>
#include <vtkPolyDataAlgorithm.h>

namespace GeoLib
{
class Point;
}

/**
 * \brief A BaseItem contains additional Information about a subtree in the StationTreeModel.
 *
 * It is used for list names in the StationTreeModel and it contains the
 * vtkObject for visualisation of the whole list in 3D.
 */
class BaseItem
{
public:
    explicit BaseItem(const QString& listName,
                      const std::vector<GeoLib::Point*>* stations = nullptr)
        : stations_(stations), vtkSource_(VtkStationSource::New())
    {
        // create the vtk-object for 3d-visualisation of this list
        static_cast<VtkStationSource*>(vtkSource_)->setStations(stations);
        static_cast<VtkStationSource*>(vtkSource_)->SetName(listName);
    }

    ~BaseItem()
    {
        vtkSource_->Delete();
    }

    /// Returns the associated QModelIndex which belongs to a Qt model
    QModelIndex modelIndex() const { return modelIndex_; }

    /// Sets the model index
    void setModelIndex( QModelIndex index ) { modelIndex_ = index; }

    const std::vector<GeoLib::Point*>* getStations() { return stations_; }

    /// Returns the Vtk polydata source object
    vtkPolyDataAlgorithm* vtkSource() const { return vtkSource_; }

private:
    QModelIndex modelIndex_;
    const std::vector<GeoLib::Point*>* stations_;

    /// The Vtk data source object. This is the starting point for a Vtk data
    /// visualization pipeline.
    vtkPolyDataAlgorithm* vtkSource_;
};
