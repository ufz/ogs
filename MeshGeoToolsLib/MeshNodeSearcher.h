/**
 * \file
 * \date Oct 24, 2013
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */
#pragma once

#include <memory>
#include <vector>

// GeoLib
#include "GeoLib/Grid.h"

// MeshGeoToolsLib
#include "MeshGeoToolsLib/SearchLength.h"
#include "MeshGeoToolsLib/SearchAllNodes.h"

// forward declaration
namespace GeoLib
{
struct GeoObject;
class Point;
class Polyline;
class Surface;
}

namespace MeshLib
{
class Mesh;
class Node;
}

namespace MeshGeoToolsLib
{
class MeshNodesOnPoint;
class MeshNodesAlongPolyline;
class MeshNodesAlongSurface;
}

namespace MeshGeoToolsLib
{

/**
 * Class for searching mesh node ids along polylines or points. This ids
 * can be used to set boundary conditions, source terms, initial conditions
 * or for outputting simulation results.
 */
class MeshNodeSearcher
{
public:
    /**
     * Constructor for objects of class MeshNodeSearcher. It calculates
     * internally used search length from the given MeshLib::Mesh.
     * @param mesh Run search algorithm on this mesh. It is assumed
     * that the mesh does not change its geometry.
     * @param search_length_algorithm Algorithm to determine the search length.
     * @param search_all_nodes switch between searching all mesh nodes and
     * searching the base nodes.
     */
    MeshNodeSearcher(MeshLib::Mesh const& mesh,
                     std::unique_ptr<MeshGeoToolsLib::SearchLength>&&
                         search_length_algorithm,
                     SearchAllNodes search_all_nodes);

    virtual ~MeshNodeSearcher();

    /**
     * Searches for the nearest mesh nodes on the given geometric object (point, polyline, surface).
     * @param geoObj a GeoLib::GeoObject where the nearest mesh node is searched for
     * @return a vector of mesh node ids
     */
    std::vector<std::size_t> getMeshNodeIDs(
        GeoLib::GeoObject const& geoObj) const;

    /**
     * Finds unique mesh nodes of each of the input points.
     *
     * \return a vector of mesh node ids.
     */
    std::vector<std::size_t> getMeshNodeIDs(
        std::vector<MathLib::Point3dWithID*> const& points) const;

    /**
     * Returns a (possibly new) mesh node searcher for the mesh.
     * A new one will be created, if it does not already exists.
     */
    static MeshNodeSearcher const& getMeshNodeSearcher(
        MeshLib::Mesh const& mesh,
        std::unique_ptr<MeshGeoToolsLib::SearchLength>&&
            search_length_algorithm);

private:
    MeshLib::Mesh const& _mesh;
    GeoLib::Grid<MeshLib::Node> _mesh_grid;
    std::unique_ptr<MeshGeoToolsLib::SearchLength> _search_length_algorithm;
    SearchAllNodes _search_all_nodes;
    // with newer compiler we can omit to use a pointer here
    mutable std::vector<MeshNodesOnPoint*> _mesh_nodes_on_points;
    mutable std::vector<MeshNodesAlongPolyline*> _mesh_nodes_along_polylines;
    mutable std::vector<MeshNodesAlongSurface*> _mesh_nodes_along_surfaces;

    /// Mesh node searcher for the meshes indexed by the meshs' ids.
    static std::vector<std::unique_ptr<MeshNodeSearcher>> _mesh_node_searchers;
};

} // end namespace MeshGeoToolsLib
