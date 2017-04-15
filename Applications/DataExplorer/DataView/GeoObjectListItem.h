/**
 * \file
 * \author Karsten Rink
 * \date   2011-02-07
 * \brief  Definition of the GeoObjectListItem class.
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 */
#pragma once

#include "TreeItem.h"

#include "GeoType.h"

#include "VtkPointsSource.h"
#include "VtkPolylinesSource.h"
#include "VtkSurfacesSource.h"
#include <QModelIndex>
#include <vtkPolyDataAlgorithm.h>

/**
 * Creates parent items for geometric object lists (i.e. points, polylines and surfaces)
 * for the GeoTreeModel
 *
 * \sa GeoTreeModel, GeoTreeItem
 */
class GeoObjectListItem : public TreeItem
{
public:
    /// Constructor for the TreeItem specifying the "Points"-subtree of a geometry.
    GeoObjectListItem(const QList<QVariant>& data,
                      TreeItem* parent,
                      const std::vector<GeoLib::Point*>* geo_data = nullptr,
                      GeoLib::GEOTYPE type = GeoLib::GEOTYPE::POINT)
        : TreeItem(data, parent),
          _vtkSource(VtkPointsSource::New()),
          _type(type)
    {
        QString geo_name = parent->data(0).toString();
        static_cast<VtkPointsSource*>(_vtkSource)->setPoints(geo_data);
        static_cast<VtkPointsSource*>(_vtkSource)->SetName(geo_name + " - Points");
    }

    /// Constructor for the TreeItem specifying the "Polylines"-subtree of a geometry.
    GeoObjectListItem(const QList<QVariant>& data,
                      TreeItem* parent,
                      const std::vector<GeoLib::Polyline*>* geo_data = nullptr,
                      GeoLib::GEOTYPE type = GeoLib::GEOTYPE::POLYLINE)
        : TreeItem(data, parent),
          _vtkSource(VtkPolylinesSource::New()),
          _type(type)
    {
        QString geo_name = parent->data(0).toString();
        static_cast<VtkPolylinesSource*>(_vtkSource)->setPolylines(geo_data);
        static_cast<VtkPolylinesSource*>(_vtkSource)->SetName(geo_name + " - Polylines");
    }

    /// Constructor for the TreeItem specifying the "Surfaces"-subtree of a geometry.
    GeoObjectListItem(const QList<QVariant>& data,
                      TreeItem* parent,
                      const std::vector<GeoLib::Surface*>* geo_data = nullptr,
                      GeoLib::GEOTYPE type = GeoLib::GEOTYPE::SURFACE)
        : TreeItem(data, parent),
          _vtkSource(VtkSurfacesSource::New()),
          _type(type)
    {
        QString geo_name = parent->data(0).toString();
        static_cast<VtkSurfacesSource*>(_vtkSource)->setSurfaces(geo_data);
        static_cast<VtkSurfacesSource*>(_vtkSource)->SetName(geo_name + " - Surfaces");
    }

    ~GeoObjectListItem() override { _vtkSource->Delete(); }
    /// Returns the type of geo-objects contained in the subtree of this item.
    GeoLib::GEOTYPE getType() { return _type; }

    /// Returns the Vtk polydata source object
    vtkPolyDataAlgorithm* vtkSource() const { return _vtkSource; }

private:
    /// The Vtk data source object. This is the starting point for a Vtk data
    /// visualization pipeline.
    vtkPolyDataAlgorithm* _vtkSource;

    GeoLib::GEOTYPE _type;
};
