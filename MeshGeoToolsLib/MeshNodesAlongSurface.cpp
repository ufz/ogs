/**
 * \author Norihiro Watanabe
 * \date   2014-03-14
 *
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MeshNodesAlongSurface.h"

#include <algorithm>

#include "BaseLib/quicksort.h"
#include "MathLib/MathTools.h"
#include "GeoLib/Surface.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"

namespace MeshGeoToolsLib
{
MeshNodesAlongSurface::MeshNodesAlongSurface(MeshLib::Mesh const& mesh,
                                             GeoLib::Surface const& sfc,
                                             double epsilon_radius,
                                             SearchAllNodes search_all_nodes)
    : mesh_(mesh), sfc_(sfc)
{
    auto& mesh_nodes = mesh_.getNodes();
    const std::size_t n_nodes(search_all_nodes == SearchAllNodes::Yes
                                  ? mesh_.getNumberOfNodes()
                                  : mesh_.getNumberOfBaseNodes());
    // loop over all nodes
    for (std::size_t i = 0; i < n_nodes; i++) {
        auto* node = mesh_nodes[i];
        if (!sfc.isPntInBoundingVolume(*node, epsilon_radius))
        {
            continue;
        }
        if (sfc.isPntInSfc(*node, epsilon_radius)) {
            msh_node_ids_.push_back(node->getID());
        }
    }
}

MeshLib::Mesh const& MeshNodesAlongSurface::getMesh () const
{
    return mesh_;
}

std::vector<std::size_t> const& MeshNodesAlongSurface::getNodeIDs () const
{
    return msh_node_ids_;
}

GeoLib::Surface const& MeshNodesAlongSurface::getSurface () const
{
    return sfc_;
}

} // end namespace MeshGeoToolsLib
