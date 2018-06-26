/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <algorithm>
#include <vector>

#include "BaseLib/makeVectorUnique.h"
#include "MeshLib/Node.h"

#include "Element.h"
#include "FaceRule.h"

namespace MeshLib
{
/// Returns a vector of node pointers containing the base nodes of the elements
/// input vector.
inline std::vector<Node*> getBaseNodes(std::vector<Element*> const& elements)
{
    std::vector<Node*> base_nodes;
    base_nodes.reserve(elements.size() * 2);  // Save some of the realloctions.

    for (auto* const e : elements)
    {
        std::copy(e->getNodes(), e->getNodes() + e->getNumberOfBaseNodes(),
                  std::back_inserter(base_nodes));
    }

    BaseLib::makeVectorUnique(base_nodes, [](Node const* a, Node* b) {
        return a->getID() < b->getID();
    });

    return base_nodes;
}

inline MathLib::Vector3 calculateNormalizedSurfaceNormal(
    MeshLib::Element const& surface_element,
    MeshLib::Element const& bulk_element)
{
    MathLib::Vector3 surface_element_normal;
    if (surface_element.getDimension() < 2)
    {
        auto const bulk_element_normal = MeshLib::FaceRule::getSurfaceNormal(
            &bulk_element);
        MathLib::Vector3 const edge_vector(*surface_element.getNode(0),
                                           *surface_element.getNode(1));
        surface_element_normal =
            MathLib::crossProduct(bulk_element_normal, edge_vector);
    }
    else
    {
        surface_element_normal =
            MeshLib::FaceRule::getSurfaceNormal(&surface_element);
    }

    surface_element_normal.normalize();
    // At the moment (2018-04-26) the surface normal is not oriented
    // according to the right hand rule
    // for correct results it is necessary to multiply the normal with
    // -1
    surface_element_normal *= -1;

    return surface_element_normal;
}

}  // namespace MeshLib
