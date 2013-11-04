/**
 * @file
 * @author git blame MeshNodesAlongPolyline.h
 * @date Aug 9, 2010
 * @brief
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef MESHNODESALONGPOLYLINE_H_
#define MESHNODESALONGPOLYLINE_H_

#include <vector>

// GeoLib
#include "Polyline.h"

// MeshLib
#include "Node.h"

namespace MeshGeoTools
{
/**
 * This class computes the ids of the mesh nodes along a polyline.
 *
 * The mesh nodes are sorted as follow:
 * [ ids of sorted nodes according to their distance to the starting point of ply ]
 */
class MeshNodesAlongPolyline
{
public:
	MeshNodesAlongPolyline(std::vector<MeshLib::Node*> const& mesh_nodes, GeoLib::Polyline const& ply, double epsilon_radius);
	std::vector<size_t> const& getNodeIDs () const;
	GeoLib::Polyline const& getPolyline () const;
	std::vector<double> const & getDistOfProjNodeFromPlyStart() const;

private:
	GeoLib::Polyline const& _ply;
	std::vector<std::size_t> _msh_node_ids;
	std::vector<double> _dist_of_proj_node_from_ply_start;
};
} // end namespace MeshGeoToolsLib

#endif /* MESHNODESALONGPOLYLINE_H_ */
