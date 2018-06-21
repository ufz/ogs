/**
 * @copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */
#pragma once

#include <vector>

namespace GeoLib
{
class Polyline;
}

namespace MeshLib
{
class Mesh;
class Element;
}

namespace MeshGeoToolsLib
{
class MeshNodeSearcher;

/**
 * This class collects element edges located along a polyline.
 * Note that internal edges are not collected in this class.
 */
class BoundaryElementsAlongPolyline
{
public:
    /**
     * Constructor
     * @param mesh             a mesh object
     * @param mshNodeSearcher  a MeshNodeSearcher object which is internally used to search mesh nodes
     * @param ply              a polyline object where edges are searched
     */
    BoundaryElementsAlongPolyline(MeshLib::Mesh const& mesh,
                                  MeshNodeSearcher const& mshNodeSearcher,
                                  GeoLib::Polyline const& ply);

    /// destructor
    virtual ~BoundaryElementsAlongPolyline();

    /// return the mesh object
    MeshLib::Mesh const& getMesh() const {return _mesh;}

    /**
     * Deploying this method the user can get access to the underlying
     * GeoLib::Polyline.
     * @return the underlying GeoLib::Polyline
     */
    GeoLib::Polyline const& getPolyline () const {return _ply;}

    /**
     * Return the vector of boundary elements (i.e. edges). The elements are sorted
     * according to their distance to the starting point of the given polyline.
     */
    std::vector<MeshLib::Element*> const& getBoundaryElements() const {return _boundary_elements; }

private:
    /**
     * Check if a vector of node IDs includes all nodes of a given element
     * @param vec_node_ids         a vector of Node IDs
     * @param edge                 Edge object whose node IDs are checked
     * @param edge_node_distances  a vector of distances of the edge nodes from the beginning of the given node ID vector
     * @return true if all element nodes are included in the vector
     */
    bool includesAllEdgeNodeIDs(const std::vector<std::size_t> &vec_node_ids, const MeshLib::Element &edge, std::vector<std::size_t> &edge_node_distances) const;

    /**
     * Modify node ordering of an edge so that its first node is closer to the beginning of a polyline than others
     * @param edge                           Element object
     * @param ply                            Polyline object
     * @param edge_node_distances_along_ply  A vector of current edge node distances along poly
     * @param node_ids_on_poly               A vector of node IDs along the polyine
     * @return A pointer to the new modified edge object. A pointer to the original edge is returned if the modification is unnecessary.
     */
    MeshLib::Element* modifyEdgeNodeOrdering(const MeshLib::Element &edge, const GeoLib::Polyline &ply, const std::vector<std::size_t> &edge_node_distances_along_ply, const std::vector<std::size_t> &node_ids_on_poly) const;

    MeshLib::Mesh const& _mesh;
    GeoLib::Polyline const& _ply;
    std::vector<MeshLib::Element*> _boundary_elements;
};

} // end namespace MeshGeoToolsLib
