/**
 * @copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "BoundaryElementsAlongPolyline.h"

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
	// search elements near the polyline
	auto node_ids_on_poly = mshNodeSearcher.getMeshNodeIDsAlongPolyline(ply);
	auto ele_ids_near_poly = MeshLib::getConnectedElementIDs(_mesh, node_ids_on_poly);

	// get a list of edges made of the nodes
	for (auto ele_id : ele_ids_near_poly) {
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
		// find edges on polyline
		for (unsigned i=0; i<e->getNEdges(); i++) {
			auto* edge = e->getEdge(i);
			// check
			size_t cnt_match = 0;
			for (size_t j=0; j<edge->getNNodes(); j++) {
				if (std::find(node_ids_on_poly.begin(), node_ids_on_poly.end(), edge->getNodeIndex(j)) != node_ids_on_poly.end())
					cnt_match++;
				else
					break;
			}
			// update the list
			if (cnt_match==edge->getNNodes())
				_boundary_elements.push_back(const_cast<MeshLib::Element*>(edge));
		}
	}
}

BoundaryElementsAlongPolyline::~BoundaryElementsAlongPolyline()
{
	for (auto p : _boundary_elements)
		delete p;
}

} // end namespace MeshGeoTools

