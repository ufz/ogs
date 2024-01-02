/**
 * \file
 * \date Aug 9, 2010
 *
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshNodesAlongPolyline.h"

#include <algorithm>

#include "BaseLib/quicksort.h"
#include "GeoLib/Polyline.h"
#include "MathLib/MathTools.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace MeshGeoToolsLib
{
MeshNodesAlongPolyline::MeshNodesAlongPolyline(MeshLib::Mesh const& mesh,
                                               GeoLib::Polyline const& ply,
                                               double epsilon_radius,
                                               SearchAllNodes search_all_nodes)
    : _mesh(mesh), _ply(ply)
{
    assert(epsilon_radius > 0);
    const std::size_t n_nodes(search_all_nodes == SearchAllNodes::Yes
                                  ? _mesh.getNumberOfNodes()
                                  : _mesh.computeNumberOfBaseNodes());
    auto& mesh_nodes = _mesh.getNodes();
    // loop over all nodes
    for (std::size_t i = 0; i < n_nodes; i++)
    {
        double dist =
            _ply.getDistanceAlongPolyline(*mesh_nodes[i], epsilon_radius);
        if (dist >= 0.0)
        {
            _msh_node_ids.push_back(mesh_nodes[i]->getID());
            _dist_of_proj_node_from_ply_start.push_back(dist);
        }
    }

    // sort the nodes along the polyline according to their distances
    BaseLib::quicksort<double>(_dist_of_proj_node_from_ply_start, 0,
                               _dist_of_proj_node_from_ply_start.size(),
                               _msh_node_ids);
}

MeshLib::Mesh const& MeshNodesAlongPolyline::getMesh() const
{
    return _mesh;
}

std::vector<std::size_t> const& MeshNodesAlongPolyline::getNodeIDs() const
{
    return _msh_node_ids;
}

GeoLib::Polyline const& MeshNodesAlongPolyline::getPolyline() const
{
    return _ply;
}
}  // end namespace MeshGeoToolsLib
