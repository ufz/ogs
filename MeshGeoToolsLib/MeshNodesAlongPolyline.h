/**
 * @file
 * @date Aug 9, 2010
 * @brief
 *
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#pragma once

#include <vector>

#include "MeshGeoToolsLib/SearchAllNodes.h"

namespace GeoLib
{
class Polyline;
}

namespace MeshLib
{
class Mesh;
}

namespace MeshGeoToolsLib
{
/**
 * This class computes the ids of the mesh nodes along a polyline.
 *
 * The mesh nodes are sorted as follow:
 * [ ids of sorted nodes according to their distance to the starting point of ply ]
 */
class MeshNodesAlongPolyline
{
public:
    /**
     * Constructor of object, that search mesh nodes along a
     * GeoLib::Polyline polyline within a given search radius. So the polyline
     * is something like a tube.
     * @param mesh Mesh the search will be performed on.
     * @param ply Along the GeoLib::Polyline ply the mesh nodes are searched.
     * @param epsilon_radius Search / tube radius
     * @param search_all_nodes switch between searching all mesh nodes and
     * searching the base nodes.
     */
    MeshNodesAlongPolyline(
        MeshLib::Mesh const& mesh, GeoLib::Polyline const& ply,
        double epsilon_radius,
        SearchAllNodes search_all_nodes);

    /// return the mesh object
    MeshLib::Mesh const& getMesh() const;

    /**
     * Access the vector of mesh node ids.
     * @return The vector of mesh node ids calculated in the constructor
     */
    std::vector<std::size_t> const& getNodeIDs () const;

    /**
     * Deploying this method the user can get access to the underlying
     * GeoLib::Polyline.
     * @return the underlying GeoLib::Polyline
     */
    GeoLib::Polyline const& getPolyline () const;

    /// return a vector of node distance from the polyline start. The vector
    /// corresponds to a vector returned in getNodeIDs()
    std::vector<double> const & getDistOfProjNodeFromPlyStart() const;

private:
    MeshLib::Mesh const& _mesh;
    GeoLib::Polyline const& _ply;
    std::vector<std::size_t> _msh_node_ids;
    std::vector<double> _dist_of_proj_node_from_ply_start;
};
} // end namespace MeshGeoToolsLib
