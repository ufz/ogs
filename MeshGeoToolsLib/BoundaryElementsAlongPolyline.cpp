/**
 * @copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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
#include "MeshLib/MeshSearch/ElementSearch.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"

namespace MeshGeoToolsLib
{

BoundaryElementsAlongPolyline::BoundaryElementsAlongPolyline(MeshLib::Mesh const& mesh, MeshNodeSearcher &mshNodeSearcher, GeoLib::Polyline const& ply)
: _mesh(mesh), _ply(ply)
{
    // search nodes and elements located along the polyline
    auto node_ids_on_poly = mshNodeSearcher.getMeshNodeIDsAlongPolyline(ply);
    MeshLib::ElementSearch es(_mesh);
    es.searchByNodeIDs(node_ids_on_poly);
    auto &ele_ids_near_ply = es.getSearchedElementIDs();

    // check all edges of the elements near the polyline
    for (auto ele_id : ele_ids_near_ply) {
        auto* e = _mesh.getElement(ele_id);
        // skip line elements
        if (e->getDimension() == 1)
            continue;
        // skip internal elements
        if (!e->isBoundaryElement())
            continue;
        // find edges on the polyline
        for (unsigned i=0; i<e->getNumberOfEdges(); i++) {
            auto* edge = e->getEdge(i);
            // check if all edge nodes are along the polyline (if yes, store a distance)
            std::vector<std::size_t> edge_node_distances_along_ply;
            if (includesAllEdgeNodeIDs(node_ids_on_poly, *edge, edge_node_distances_along_ply)) {
                auto* new_edge = modifyEdgeNodeOrdering(*edge, ply, edge_node_distances_along_ply, node_ids_on_poly);
                if (edge != new_edge)
                    delete edge;
                _boundary_elements.push_back(new_edge);
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

bool BoundaryElementsAlongPolyline::includesAllEdgeNodeIDs(const std::vector<std::size_t> &vec_node_ids, const MeshLib::Element &edge, std::vector<std::size_t> &edge_node_distances) const
{
    unsigned j=0;
    for (; j<edge.getNumberOfBaseNodes(); j++) {
        auto itr = std::find(vec_node_ids.begin(), vec_node_ids.end(), edge.getNodeIndex(j));
        if (itr != vec_node_ids.end())
            edge_node_distances.push_back(std::distance(vec_node_ids.begin(), itr));
        else
            break;
    }
    return (j==edge.getNumberOfBaseNodes());
}

MeshLib::Element* BoundaryElementsAlongPolyline::modifyEdgeNodeOrdering(const MeshLib::Element &edge, const GeoLib::Polyline &ply, const std::vector<std::size_t> &edge_node_distances_along_ply, const std::vector<std::size_t> &node_ids_on_poly) const
{
    MeshLib::Element* new_edge = const_cast<MeshLib::Element*>(&edge);
    // The first node of the edge should be always closer to the beginning of the polyline than other nodes.
    // Otherwise, create a new element with reversed local node index
    if (edge_node_distances_along_ply.front() > edge_node_distances_along_ply.back()
            || (ply.isClosed() && edge_node_distances_along_ply.back() == node_ids_on_poly.size()-1)) {
        MeshLib::Node** new_nodes = new MeshLib::Node*[edge.getNumberOfBaseNodes()];
        std::reverse_copy(edge.getNodes(), edge.getNodes()+edge.getNumberOfBaseNodes(), new_nodes);
        new_edge = new MeshLib::Line(new_nodes);
    }
    return new_edge;
}

} // end namespace MeshGeoToolsLib

