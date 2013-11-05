/**
 * @file
 * @date Oct 24, 2013
 * @brief
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */
#ifndef MESHNODESEARCHER_H_
#define MESHNODESEARCHER_H_

#include <vector>

// GeoLib
#include "Point.h"
#include "Polyline.h"
#include "Grid.h"

// MeshLib
#include "Mesh.h"
#include "Node.h"

// forward declaration
namespace MeshGeoTools
{
class MeshNodesAlongPolyline;
}

namespace MeshGeoTools
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
	MeshNodeSearcher(MeshLib::Mesh const& mesh);
	virtual ~MeshNodeSearcher();

	/**
	 * Searchs for the node nearest by the given point. If there are two nodes
	 * with the same distance the id of the one that was first found will be
	 * returned. The algorithm for the search is using GeoLib::Grid data
	 * structure.
	 * @param pnt a GeoLib::Point the nearest mesh node is searched for
	 * @return the id of the nearest mesh node
	 */
	std::size_t getMeshNodeIDForPoint(GeoLib::Point const& pnt) const;

	/**
	 * Searchs for the nearest mesh nodes along a GeoLib::Polyline.
	 * The search for mesh nodes along a specific polyline will be performed
	 * only once. The result ids will be stored inside an object
	 * (@see class MeshGeoTools::MeshNodesAlongPolyline).
	 * @param ply the GeoLib::Polyline the nearest mesh nodes are searched for
	 * @return a vector of mesh node ids
	 */
	std::vector<std::size_t> const& getMeshNodeIDsAlongPolyline(GeoLib::Polyline const& ply);

private:
	MeshLib::Mesh const& _mesh;
	GeoLib::Grid<MeshLib::Node> _mesh_grid;
	double _search_length;
	// with newer compiler we can omit to use a pointer here
	std::vector<MeshNodesAlongPolyline*> _mesh_nodes_along_polylines;
};

} // end namespace MeshGeoTools

#endif /* MESHNODESEARCHER_H_ */
