/**
 * \file
 * \author Karsten Rink
 * \date   2015-04-10
 * \brief  Definition of the projectMeshOntoPlane
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "MathLib/MathTools.h"
#include "MathLib/Point3d.h"
#include "MathLib/Vector3.h"

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/MeshEditing/DuplicateMeshComponents.h"


namespace MeshLib {

/**
 * Projects all nodes of a mesh onto a plane specified by a point of origin and a normal vector.
 * Overlapping elements, collapsed nodes, and other issues are not handled by the method.
 * The normal vector need not be normalized.
 */
inline
MeshLib::Mesh* projectMeshOntoPlane(MeshLib::Mesh const& mesh,
                                    MathLib::Point3d const& plane_origin,
                                    MathLib::Vector3 const& plane_normal)
{
    std::size_t const n_nodes (mesh.getNumberOfNodes());
    std::vector<MeshLib::Node*> const& nodes (mesh.getNodes());
    Eigen::Vector3d normal({plane_normal[0], plane_normal[1], plane_normal[2]});
    normal.normalize();
    std::vector<MeshLib::Node*> new_nodes;
    new_nodes.reserve(n_nodes);
    for (std::size_t i=0; i<n_nodes; ++i)
    {
        auto const node =
            Eigen::Map<Eigen::Vector3d const>(nodes[i]->getCoords());
        Eigen::Vector3d const v =
            node - Eigen::Map<Eigen::Vector3d const>(plane_origin.getCoords());
        double const dist(v.dot(normal));
        Eigen::Vector3d const new_node = node - dist * normal;
        new_nodes.push_back(
            new MeshLib::Node(new_node[0], new_node[1], new_node[2]));
    }

    return new MeshLib::Mesh("Projected_Mesh", new_nodes,
                             MeshLib::copyElementVector(mesh.getElements(), new_nodes));
}

} // end namespace MeshLib
