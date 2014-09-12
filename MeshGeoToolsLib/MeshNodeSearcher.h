/**
 * @file
 * @date Oct 24, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */
#ifndef MESHNODESEARCHER_H_
#define MESHNODESEARCHER_H_

#include <vector>

#include "boost/optional.hpp"

// GeoLib
#include "Point.h"
#include "Polyline.h"
#include "Grid.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"

// forward declaration
namespace MeshGeoToolsLib
{
class MeshNodesAlongPolyline;
class MeshNodesAlongSurface;
}

namespace MeshGeoToolsLib
{
/**
 * Class for searching mesh node ids along polylines or points. This ids
 * can be used to set boundary conditions, source terms, initial conditions
 * or for outputting simulation results.
 */
class MeshNodeSearcher
{
public:
	/**
	 * Constructor for objects of class MeshNodeSearcher. It calculates
	 * internal used search length out of the given MeshLib::Mesh.
	 * @param mesh The mesh within the search will be performed. It is asumed
	 * that the mesh does not change its geometry.
	 */
	explicit MeshNodeSearcher(MeshLib::Mesh const& mesh);
	virtual ~MeshNodeSearcher();

	/**
	 * Searches for the nearest mesh nodes on the given geometric object (point, polyline, surface).
	 * @param geoObj a GeoLib::GeoObject where the nearest mesh node is searched for
	 * @return a vector of mesh node ids
	 */
	std::vector<std::size_t> getMeshNodeIDs(GeoLib::GeoObject const& geoObj);

	/**
	 * Searches for the node nearest by the given point. If there are two nodes
	 * with the same distance the id of the one that was first found will be
	 * returned. The algorithm for the search is using GeoLib::Grid data
	 * structure.
	 * @param pnt a GeoLib::Point the nearest mesh node is searched for
	 * @return the id of the nearest mesh node
	 */
	boost::optional<std::size_t> getMeshNodeIDForPoint(GeoLib::Point const& pnt) const;

	/**
	 * Searches for the nearest mesh nodes along a GeoLib::Polyline.
	 * The search for mesh nodes along a specific polyline will be performed
	 * only once. The result ids will be stored inside an object
	 * (@see class MeshGeoTools::MeshNodesAlongPolyline).
	 * @param ply the GeoLib::Polyline the nearest mesh nodes are searched for
	 * @return a vector of mesh node ids
	 */
	std::vector<std::size_t> const& getMeshNodeIDsAlongPolyline(GeoLib::Polyline const& ply);

	/**
	 * Searches for the nearest mesh nodes along a GeoLib::Surface.
	 * The search for mesh nodes along a specific surface will be performed
	 * only once. The result ids will be stored inside an object
	 * (@see class MeshGeoTools::MeshNodesAlongSurface).
	 * @param sfc the GeoLib::Surface the nearest mesh nodes are searched for
	 * @return a vector of mesh node ids
	 */
	std::vector<std::size_t> const& getMeshNodeIDsAlongSurface(GeoLib::Surface const& sfc);

	/**
	 * Return a MeshNodesAlongPolyline object for the given GeoLib::Polyline object.
	 * @param ply the GeoLib::Polyline the nearest mesh nodes are searched for
	 * @return a reference to a MeshNodesAlongPolyline object
	 */
	MeshNodesAlongPolyline& getMeshNodesAlongPolyline(GeoLib::Polyline const& ply);

	/**
	 * Return a MeshNodesAlongSurface object for the given GeoLib::Surface object.
	 * @param sfc the GeoLib::Surface the nearest mesh nodes are searched for
	 * @return a reference to a MeshNodesAlongSurface object
	 */
	MeshNodesAlongSurface& getMeshNodesAlongSurface(GeoLib::Surface const& sfc);

private:
	MeshLib::Mesh const& _mesh;
	GeoLib::Grid<MeshLib::Node> _mesh_grid;
	double _search_length;
	// with newer compiler we can omit to use a pointer here
	std::vector<MeshNodesAlongPolyline*> _mesh_nodes_along_polylines;
	std::vector<MeshNodesAlongSurface*> _mesh_nodes_along_surfaces;
};

} // end namespace MeshGeoTools

#endif /* MESHNODESEARCHER_H_ */
