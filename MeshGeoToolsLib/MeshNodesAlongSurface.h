/**
 * \author Norihiro Watanabe
 * \date   2014-03-14
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "MeshGeoToolsLib/SearchAllNodes.h"

namespace GeoLib
{
class Surface;
}

namespace MeshLib
{
class Mesh;
}

namespace MeshGeoToolsLib
{
/**
 * This class computes the ids of the mesh nodes along a surface.
 */
class MeshNodesAlongSurface
{
public:
    /**
     * Constructor of object, that search mesh nodes along a
     * GeoLib::Surface object within a given search radius.
     * @param mesh Mesh the search will be performed on.
     * @param sfc Along the GeoLib::Surface sfc the mesh nodes are searched.
     * @param epsilon Euclidean distance tolerance value. Is the distance
     * between a mesh node and the surface smaller than that value it is a mesh
     * node near the surface.
     * @param search_all_nodes switch between searching all mesh nodes and
     * searching the base nodes.
     */
    MeshNodesAlongSurface(MeshLib::Mesh const& mesh, GeoLib::Surface const& sfc,
                          double epsilon, SearchAllNodes search_all_nodes);

    /// return the mesh object
    MeshLib::Mesh const& getMesh() const;

    /**
     * Access the vector of mesh node ids.
     * @return The vector of mesh node ids calculated in the constructor
     */
    std::vector<std::size_t> const& getNodeIDs () const;

    /**
     * Deploying this method the user can get access to the underlying
     * GeoLib::Surface.
     * @return the underlying GeoLib::Surface
     */
    GeoLib::Surface const& getSurface () const;

private:
    MeshLib::Mesh const& _mesh;
    GeoLib::Surface const& _sfc;
    std::vector<std::size_t> _msh_node_ids;
};
} // end namespace MeshGeoToolsLib
