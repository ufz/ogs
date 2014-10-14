/**
 * \author Norihiro Watanabe
 * \date   2014-03-14
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESHNODESALONGSURFACE_H_
#define MESHNODESALONGSURFACE_H_

#include <vector>

// GeoLib
#include "Surface.h"

// MeshLib
#include "Node.h"

namespace MeshGeoToolsLib
{
/**
 * This class computes the ids of the mesh nodes along a surface.
 */
class MeshNodesAlongSurface
{
public:
	/**
	 * Constructor of object, that search mesh nodes along a
	 * GeoLib::Surface object within a given search radius.
	 * @param mesh_nodes Nodes the search will be performed on.
	 * @param sfc Along the GeoLib::Surface sfc the mesh nodes are searched.
	 */
	MeshNodesAlongSurface(MeshLib::Mesh const& mesh, GeoLib::Surface const& sfc, bool search_all_nodes = true);

	/// return the mesh object
	MeshLib::Mesh const& getMesh() const;

	/**
	 * Access the vector of mesh node ids.
	 * @return The vector of mesh node ids calculated in the constructor
	 */
	std::vector<std::size_t> const& getNodeIDs () const;

	/**
	 * Deploying this method the user can get access to the underlying
	 * GeoLib::Surface.
	 * @return the underlying GeoLib::Surface
	 */
	GeoLib::Surface const& getSurface () const;

private:
	MeshLib::Mesh const& _mesh;
	GeoLib::Surface const& _sfc;
	std::vector<std::size_t> _msh_node_ids;
};
} // end namespace MeshGeoToolsLib

#endif /* MESHNODESALONGSURFACE_H_ */
