/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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
    MeshLib::Mesh const &mesh, MeshNodeSearcher &mshNodeSearcher,
    GeoLib::Point const &point)
    : _mesh(mesh), _point(point)
{
    auto const node_ids = mshNodeSearcher.getMeshNodeIDs(_point);
    assert(node_ids.size() == 1);
    std::array<MeshLib::Node*, 1> const nodes = {{
        const_cast<MeshLib::Node*>(_mesh.getNode(node_ids[0]))}};

    _boundary_elements.push_back(new MeshLib::Point{nodes, node_ids[0]});
}

BoundaryElementsAtPoint::~BoundaryElementsAtPoint()
{
    for (auto p : _boundary_elements)
        delete p;
}

}  // MeshGeoToolsLib
