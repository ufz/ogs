/**
 * \file
 * \date Aug 9, 2010
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshNodesAlongPolyline.h"

#include <algorithm>

#include "BaseLib/quicksort.h"
#include "MathLib/MathTools.h"
#include "GeoLib/Polyline.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace MeshGeoToolsLib
{
MeshNodesAlongPolyline::MeshNodesAlongPolyline(MeshLib::Mesh const& mesh,
                                               GeoLib::Polyline const& ply,
                                               double epsilon_radius,
                                               SearchAllNodes search_all_nodes)
    : mesh_(mesh), ply_(ply)
{
    assert(epsilon_radius > 0);
    const std::size_t n_nodes(search_all_nodes == SearchAllNodes::Yes
                                  ? mesh_.getNumberOfNodes()
                                  : mesh_.getNumberOfBaseNodes());
    auto &mesh_nodes = mesh_.getNodes();
    // loop over all nodes
    for (std::size_t i = 0; i < n_nodes; i++) {
        double dist = ply_.getDistanceAlongPolyline(*mesh_nodes[i], epsilon_radius);
        if (dist >= 0.0) {
            msh_node_ids_.push_back(mesh_nodes[i]->getID());
            dist_of_proj_node_from_ply_start_.push_back(dist);
        }
    }

    // sort the nodes along the polyline according to their distances
    BaseLib::quicksort<double> (dist_of_proj_node_from_ply_start_, 0,
                    dist_of_proj_node_from_ply_start_.size(), msh_node_ids_);
}

MeshLib::Mesh const& MeshNodesAlongPolyline::getMesh () const
{
    return mesh_;
}

std::vector<std::size_t> const& MeshNodesAlongPolyline::getNodeIDs () const
{
    return msh_node_ids_;
}

GeoLib::Polyline const& MeshNodesAlongPolyline::getPolyline () const
{
    return ply_;
}
}  // end namespace MeshGeoToolsLib
