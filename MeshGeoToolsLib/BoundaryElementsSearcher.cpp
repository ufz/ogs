/**
 * @copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "BoundaryElementsSearcher.h"

#include "MeshLib/Elements/Element.h"
#include "MeshLib/MeshSearcher.h"

namespace MeshGeoToolsLib
{

BoundaryElementsSearcher::BoundaryElementsSearcher(MeshLib::Mesh const& mesh, MeshNodeSearcher &mshNodeSearcher) : _mesh(mesh), _mshNodeSearcher(mshNodeSearcher)
{}

std::vector<MeshLib::Element*> BoundaryElementsSearcher::getBoundaryElements(GeoLib::GeoObject const& geoObj)
{
	std::vector<MeshLib::Element*> vec_elements;
	switch (geoObj.getGeoType()) {
	case GeoLib::GEOTYPE::POLYLINE:
		vec_elements = this->getBoundaryElementsAlongPolyline(*dynamic_cast<const GeoLib::Polyline*>(&geoObj));
		break;
	default:
		break;
	}
	return vec_elements;
}

std::vector<MeshLib::Element*> BoundaryElementsSearcher::getBoundaryElementsAlongPolyline(GeoLib::Polyline const& ply)
{
	// serach elements near the polyline
	auto node_ids_on_poly = _mshNodeSearcher.getMeshNodeIDsAlongPolyline(ply);
	auto ele_ids_near_poly = MeshLib::getConnectedElementIDs(_mesh, node_ids_on_poly);

	// get a list of edges made of the nodes
	std::vector<MeshLib::Element*> vec_edges_on_poly;
	for (auto ele_id : ele_ids_near_poly) {
		auto* e = _mesh.getElement(ele_id);
		for (unsigned i=0; i<e->getNEdges(); i++) {
			auto* edge = e->getEdge(i);
			//TODO where should we store and delete this new object?
			//TODO avoid duplicated entries
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
				vec_edges_on_poly.push_back(const_cast<MeshLib::Element*>(edge));
		}
	}
	return vec_edges_on_poly;
}


} // end namespace MeshGeoTools

