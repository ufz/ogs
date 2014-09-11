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
#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshSearcher.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"

namespace MeshGeoToolsLib
{

BoundaryElementsAlongPolyline::BoundaryElementsAlongPolyline(MeshLib::Mesh const& mesh, MeshNodeSearcher &mshNodeSearcher, GeoLib::Polyline const& ply)
: _mesh(mesh), _ply(ply)
{
	// search nodes located along the polyline
	auto node_ids_on_poly = mshNodeSearcher.getMeshNodeIDsAlongPolyline(ply);
	std::vector<std::size_t> work_node_ids(node_ids_on_poly);
	if (ply.isClosed())
		work_node_ids.push_back(work_node_ids[0]);

	// search edges
	for (unsigned i=0; i<work_node_ids.size()-1; i++) {
		std::vector<std::size_t> edge_nodeIDs = {work_node_ids[i], work_node_ids[i+1]};
		// find a shared element
		auto* node1 = mesh.getNode(edge_nodeIDs[0]);
		auto* node2 = mesh.getNode(edge_nodeIDs[1]);
		auto it = std::find_first_of(node1->getElements().begin(), node1->getElements().end(), node2->getElements().begin(), node2->getElements().end());
		assert (it!=node1->getElements().end());
		auto* e = *it;
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
		// find edges on polyline
		for (unsigned i=0; i<e->getNEdges(); i++) {
			auto* edge = e->getEdge(i);
			// check
			size_t cnt_match = 0;
			for (size_t j=0; j<edge->getNNodes(); j++) {
				if (std::find(edge_nodeIDs.begin(), edge_nodeIDs.end(), edge->getNodeIndex(j)) != edge_nodeIDs.end())
					cnt_match++;
				else
					break;
			}
			// update the list
			if (cnt_match==edge->getNNodes()) {
				_boundary_elements.push_back(const_cast<MeshLib::Element*>(edge));
				break;
			}
		}
	}
}

BoundaryElementsAlongPolyline::~BoundaryElementsAlongPolyline()
{
	for (auto p : _boundary_elements)
		delete p;
}

} // end namespace MeshGeoTools

