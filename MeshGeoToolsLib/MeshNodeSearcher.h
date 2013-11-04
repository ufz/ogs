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

class MeshNodeSearcher
{
public:
	MeshNodeSearcher(MeshLib::Mesh const& mesh);
	virtual ~MeshNodeSearcher();

	std::size_t getMeshNodeIDForPoint(GeoLib::Point const& pnt) const;
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
