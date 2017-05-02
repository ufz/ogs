/**
 * \file
 * \author Karsten Rink
 * \date   2010-01-20
 * \brief  Definition of the BaseItem class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#include "Point.h"

#include "VtkStationSource.h"
#include <QModelIndex>
#include <vtkPolyDataAlgorithm.h>

/**
 * \brief A BaseItem contains additional Information about a subtree in the StationTreeModel.
 *
 * It is used for list names in the StationTreeModel and it contains the
 * vtkObject for visualisation of the whole list in 3D.
 */
class BaseItem
{
public:
    BaseItem(const QString& listName,
             const std::vector<GeoLib::Point*>* stations = nullptr)
        : _stations(stations), _vtkSource(VtkStationSource::New())
    {
        // create the vtk-object for 3d-visualisation of this list
        static_cast<VtkStationSource*>(_vtkSource)->setStations(stations);
        static_cast<VtkStationSource*>(_vtkSource)->SetName(listName);
    }

    ~BaseItem()
    {
        _vtkSource->Delete();
    }

    /// Returns the associated QModelIndex which belongs to a Qt model
    QModelIndex modelIndex() const { return _modelIndex; }

    /// Sets the model index
    void setModelIndex( QModelIndex index ) { _modelIndex = index; }

    const std::vector<GeoLib::Point*>* getStations() { return _stations; }

    /// Returns the Vtk polydata source object
    vtkPolyDataAlgorithm* vtkSource() const { return _vtkSource; }

private:
    QModelIndex _modelIndex;
    const std::vector<GeoLib::Point*>* _stations;

    /// The Vtk data source object. This is the starting point for a Vtk data
    /// visualization pipeline.
    vtkPolyDataAlgorithm* _vtkSource;
};
