/**
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "BoundaryElementsAtPoint.h"

#include "GeoLib/Point.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Elements/Point.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshSearch/ElementSearch.h"
#include "MeshLib/Node.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"

namespace MeshGeoToolsLib
{
BoundaryElementsAtPoint::BoundaryElementsAtPoint(
    MeshLib::Mesh const& mesh, MeshNodeSearcher const& mshNodeSearcher,
    GeoLib::Point const& point)
    : _mesh(mesh), _point(point)
{
    auto const node_ids = mshNodeSearcher.getMeshNodeIDs(_point);
    if (node_ids.empty())
        OGS_FATAL(
            "BoundaryElementsAtPoint: the mesh node searcher was unable to "
            "locate the point (%f, %f, %f) in the mesh.",
            _point[0], _point[1], _point[2]);
    if (node_ids.size() > 1)
        OGS_FATAL(
            "BoundaryElementsAtPoint: the mesh node searcher found %d points "
            "near the requested point (%f, %f, %f) in the mesh, while exactly "
            "one is expected.",
            node_ids.size(), _point[0], _point[1], _point[2]);

    auto& mesh_nodes =
        const_cast<std::vector<MeshLib::Node*>&>(_mesh.getNodes());

    std::array<MeshLib::Node*, 1> const nodes = {{
        const_cast<MeshLib::Node*>(_mesh.getNode(node_ids[0]))}};

    _boundary_elements.push_back(new MeshLib::Point{nodes, node_ids[0]});
    for (auto const* bulk_element : _mesh.getNode(node_ids[0])->getElements())
    {
        _bulk_ids.emplace_back(
            bulk_element->getID(),
            bulk_element->identifyFace(mesh_nodes.data() + node_ids[0]));
    }
}

BoundaryElementsAtPoint::~BoundaryElementsAtPoint()
{
    for (auto p : _boundary_elements)
        delete p;
}

std::vector<std::pair<std::size_t, unsigned>> const&
BoundaryElementsAtPoint::getBulkIDs() const
{
    return _bulk_ids;
}

}  // MeshGeoToolsLib
