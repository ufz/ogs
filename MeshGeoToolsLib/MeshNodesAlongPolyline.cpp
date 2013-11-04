/**
 * @file
 * @date Aug 9, 2010
 * @brief
 *
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

// BaseLib
#include "quicksort.h"

// MeshGeoToolsLib
#include "MeshNodesAlongPolyline.h"

#include "MathTools.h"

#include <algorithm>

namespace MeshGeoTools
{
MeshNodesAlongPolyline::MeshNodesAlongPolyline(
		std::vector<MeshLib::Node*> const& mesh_nodes,
		GeoLib::Polyline const& ply,
		double epsilon_radius) :
	_ply(ply)
{
	const std::size_t n_nodes (mesh_nodes.size());
	// loop over all line segments of the polyline
	for (size_t k = 0; k < _ply.getNumberOfPoints() - 1; k++) {
		double act_length_of_ply(_ply.getLength(k));
		double seg_length (sqrt(MathLib::sqrDist(_ply.getPoint(k), _ply.getPoint(k + 1))));
		double lower_lambda (- epsilon_radius / seg_length);
		double upper_lambda (1 + epsilon_radius / seg_length);

		// loop over all nodes
		for (size_t j = 0; j < n_nodes; j++) {
			double dist, lambda;

			// is the orthogonal projection of the j-th node to the
			// line g(lambda) = _ply->getPoint(k) + lambda * (_ply->getPoint(k+1) - _ply->getPoint(k))
			// at the k-th line segment of the polyline, i.e. 0 <= lambda <= 1?
			if (MathLib::calcProjPntToLineAndDists(mesh_nodes[j]->getCoords(),
							(_ply.getPoint(k))->getCoords(), (_ply.getPoint(k + 1))->getCoords(),
							lambda, dist) <= epsilon_radius) {
				if (lower_lambda <= lambda && lambda <= upper_lambda) {
					// check if node id is already in the vector
					if (std::find(_msh_node_ids.begin(), _msh_node_ids.end(),
									mesh_nodes[j]->getID()) == _msh_node_ids.end()) {
						_msh_node_ids.push_back(mesh_nodes[j]->getID());
						_dist_of_proj_node_from_ply_start.push_back(act_length_of_ply + dist);
					}
				} // end if lambda
			}
		} // end node loop
	} // end line segment loop

	// sort the nodes along the polyline according to their distances
	BaseLib::Quicksort<double> (_dist_of_proj_node_from_ply_start, 0,
					_dist_of_proj_node_from_ply_start.size(), _msh_node_ids);
}

std::vector<size_t> const& MeshNodesAlongPolyline::getNodeIDs () const
{
	return _msh_node_ids;
}

GeoLib::Polyline const& MeshNodesAlongPolyline::getPolyline () const
{
	return _ply;
}

std::vector<double> const& MeshNodesAlongPolyline::getDistOfProjNodeFromPlyStart() const
{
	return _dist_of_proj_node_from_ply_start;
}
} // end namespace MeshGeoTools
