/**
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef MESHNODESONPOINT_H_
#define MESHNODESONPOINT_H_

#include <vector>

#include "GeoLib/Point.h"
#include "GeoLib/Grid.h"

#include "MeshLib/Node.h"

namespace MeshLib
{
class Mesh;
}

namespace MeshGeoToolsLib
{
/**
 * This class computes the ids of the mesh nodes located at a given point.
 */
class MeshNodesOnPoint
{
public:
    /**
     * Constructor of object, that search mesh nodes at a GeoLib::Point point
     * within a given search radius.
     * @param mesh Mesh object whose nodes are searched
     * @param mesh_grid Grid object constructed with mesh nodes
     * @param pnt a point
     * @param epsilon_radius Search radius
     * @param search_all_nodes whether this shearches all nodes or only base nodes
     */
    MeshNodesOnPoint(MeshLib::Mesh const& mesh, GeoLib::Grid<MeshLib::Node> const &mesh_grid,
            GeoLib::Point const& pnt, double epsilon_radius, bool search_all_nodes = true);

    /// return the mesh object
    MeshLib::Mesh const& getMesh() const { return _mesh; }

    /**
     * Access the vector of mesh node ids.
     * @return The vector of mesh node ids calculated in the constructor
     */
    std::vector<std::size_t> const& getNodeIDs () const { return _msh_node_ids; }

    /**
     * Deploying this method the user can get access to the underlying
     * GeoLib::Point.
     * @return the underlying GeoLib::Point
     */
    GeoLib::Point const& getPoint () const { return _pnt; }

private:
    MeshLib::Mesh const& _mesh;
    GeoLib::Point const& _pnt;
    std::vector<std::size_t> _msh_node_ids;
};
} // end namespace MeshGeoToolsLib

#endif /* MESHNODESONPOINT_H_ */
