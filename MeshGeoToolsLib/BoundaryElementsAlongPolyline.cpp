/**
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "BoundaryElementsAlongPolyline.h"

#include <algorithm>

#include "GeoLib/Polyline.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Line.h"
#include "MeshLib/MeshSearcher.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"

namespace MeshGeoToolsLib
{

BoundaryElementsAlongPolyline::BoundaryElementsAlongPolyline(MeshLib::Mesh const& mesh, MeshNodeSearcher &mshNodeSearcher, GeoLib::Polyline const& ply)
: _mesh(mesh), _ply(ply)
{
	// search nodes and elements located along the polyline
	auto node_ids_on_poly = mshNodeSearcher.getMeshNodeIDsAlongPolyline(ply);
	auto ele_ids_near_ply = MeshLib::getConnectedElementIDs(_mesh, node_ids_on_poly);

	// check all edges of the elements near the polyline
	for (auto ele_id : ele_ids_near_ply) {
		auto* e = _mesh.getElement(ele_id);
		// skip internal elements
		bool isOuterElement = false;
		for (unsigned i=0; i<e->getNNeighbors(); i++) {
			if (!e->getNeighbor(i)) {
				isOuterElement = true;
				break;
			}
		}
		if (!isOuterElement)
			continue;
		// find edges on the polyline
		for (unsigned i=0; i<e->getNEdges(); i++) {
			auto* edge = e->getEdge(i);
			// check if an edge node is on the polyline (if yes, store a distance)
			std::vector<std::size_t> vec_matched_node_distance_along_ply;
			for (size_t j=0; j<edge->getNNodes(); j++) {
				auto itr = std::find(node_ids_on_poly.begin(), node_ids_on_poly.end(), edge->getNodeIndex(j));
				if (itr != node_ids_on_poly.end())
					vec_matched_node_distance_along_ply.push_back(std::distance(node_ids_on_poly.begin(), itr));
				else
					break;
			}
			// the edge is picked if its all nodes are on the polyline
			if (vec_matched_node_distance_along_ply.size()==edge->getNNodes()) {
				MeshLib::Element* picked_ele = const_cast<MeshLib::Element*>(edge);
				// The first node of the edge should be always closer to the beginning of the polyline than other nodes.
				// Otherwise, create a new element with reversed local node index
				if (vec_matched_node_distance_along_ply.front() > vec_matched_node_distance_along_ply.back()
						|| (ply.isClosed() && vec_matched_node_distance_along_ply.back() == node_ids_on_poly.size()-1)) {
					MeshLib::Node** new_nodes = new MeshLib::Node*[edge->getNNodes()];
					std::reverse_copy(edge->getNodes(), edge->getNodes()+edge->getNNodes(), new_nodes);
					picked_ele = new MeshLib::Line(new_nodes);
					delete edge;
				}
				_boundary_elements.push_back(picked_ele);
			} else {
				delete edge;
			}
		}
	}

	// sort picked edges according to a distance of their first node along the polyline
	std::sort(_boundary_elements.begin(), _boundary_elements.end(),
			[&](MeshLib::Element*e1, MeshLib::Element*e2)
			{
				std::size_t dist1 = std::distance(node_ids_on_poly.begin(),
					std::find(node_ids_on_poly.begin(), node_ids_on_poly.end(), e1->getNodeIndex(0)));
				std::size_t dist2 = std::distance(node_ids_on_poly.begin(),
					std::find(node_ids_on_poly.begin(), node_ids_on_poly.end(), e2->getNodeIndex(0)));
				return (dist1 < dist2);
			});
}

BoundaryElementsAlongPolyline::~BoundaryElementsAlongPolyline()
{
	for (auto p : _boundary_elements)
		delete p;
}

} // end namespace MeshGeoTools

