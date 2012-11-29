/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * \file GeoTreeModel.h
 *
 * Created on 2011-02-07 by Karsten Rink
 */

#ifndef GEOTREEMODEL_H
#define GEOTREEMODEL_H

#include <vector>

#include "GeoType.h"
#include "PointVec.h"
#include "PolylineVec.h"
#include "SurfaceVec.h"
#include "TreeModel.h"
#include <vtkPolyDataAlgorithm.h>

namespace GeoLib
{
class GeoObject;
}

class QString;
class QModelIndex;
class GeoTreeItem;
class GeoObjectListItem;

/**
 * \brief A model for the GeoTreeView implementing a tree as a double-linked list.
 * \sa TreeModel, GeoTreeView, TreeItem, GeoTreeItem
 */
class GeoTreeModel : public TreeModel
{
	Q_OBJECT

public:
	GeoTreeModel( QObject* parent = 0 );
	~GeoTreeModel();

	/**
	 * Inserts a new subtree under _rootItem for a geometry with the given name along with a subtree named "Points" for that new geometry.
	 * \param geoName Name of the new subtree. If no name is given a default name is assigned.
	 * \param pointVec The list of points to be added as children of that subtree (no geometry can be added with a set of points!)
	 */
	void addPointList(QString geoName, const GeoLib::PointVec* pointVec);
	/// Adds a subtree "Polylines" to an existing geometry with the given name.
	void addPolylineList(QString geoName, const GeoLib::PolylineVec* polylineVec);
	/// Appends polylines to the "Polyline"-subtree
	void appendPolylines(const std::string &name, const GeoLib::PolylineVec* polylineVec);
	/// Adds a subtree "Surfaces" to an existing geometry with the given name.
	void addSurfaceList(QString geoName, const GeoLib::SurfaceVec* surfaceVec);
	/// Appends surfaces to the "Surface"-subtree
	void appendSurfaces(const std::string &name, GeoLib::SurfaceVec* surfaceVec);
	/// Returns a list of all existing geometries.
	const std::vector<GeoTreeItem*> &getLists() { return _lists; }
	/**
	 * Removes a geometry (points, polylines & surfaces) or just a specified subtree indicated by type.
	 * Note that points cannot be deleted as long as other objects exist that depend on them.
	 */
	void removeGeoList(const std::string &name, GeoLib::GEOTYPE type = GeoLib::INVALID);

	void setNameForItem(const std::string &name, GeoLib::GEOTYPE type, std::size_t id, std::string item_name);

	/*
	 * Returns the geo-object specified by the given index.
	 * \param index Index of the requested item
	 * \param listName Here, the method will put the name of the geometry this object belongs to.
	 * \return A geo-object (Point / Polyline / Surface)
	 */
	//const GeoLib::GeoObject* objectFromIndex( const QModelIndex& index, QString &geoName ) const;

	/// Returns the vtk-object indicated by type of the geometry indicated by name.
	vtkPolyDataAlgorithm* vtkSource(const std::string &name, GeoLib::GEOTYPE type) const;

private:
	/// Adds children to the "Polylines" node
	void addChildren(GeoObjectListItem* plyList,
	                 const GeoLib::PolylineVec* polyline_vec,
	                 std::size_t start_index,
	                 std::size_t end_index);

	/// Adds children to the "Surfaces" node
	void addChildren(GeoObjectListItem* sfcList,
	                 const GeoLib::SurfaceVec* surface_vec,
	                 std::size_t start_index,
	                 std::size_t end_index);

	std::vector<GeoTreeItem*> _lists;
};

#endif //GEOTREEMODEL_H
