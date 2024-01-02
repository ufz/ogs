/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Geometry>
#include <algorithm>
#include <vector>

#include "BaseLib/Algorithm.h"
#include "Element.h"
#include "FaceRule.h"
#include "MeshLib/Node.h"

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

    BaseLib::makeVectorUnique(base_nodes, MeshLib::idsComparator<Node*>);

    return base_nodes;
}

inline Eigen::Vector3d calculateNormalizedSurfaceNormal(
    MeshLib::Element const& surface_element,
    MeshLib::Element const& bulk_element)
{
    Eigen::Vector3d surface_element_normal;

    switch (surface_element.getDimension())
    {
        case 2:
            surface_element_normal =
                MeshLib::FaceRule::getSurfaceNormal(surface_element);
            break;
        case 1:
        {
            auto const bulk_element_normal =
                MeshLib::FaceRule::getSurfaceNormal(bulk_element);
            auto const& v0 = surface_element.getNode(0)->asEigenVector3d();
            auto const& v1 = surface_element.getNode(1)->asEigenVector3d();
            Eigen::Vector3d const edge_vector = v1 - v0;
            surface_element_normal = -bulk_element_normal.cross(edge_vector);
            break;
        }
        case 0:
        {
            assert(surface_element.getCellType() == CellType::POINT1);
            assert(bulk_element.getCellType() == CellType::LINE2 ||
                   bulk_element.getCellType() == CellType::LINE3);

            auto const& x = surface_element.getNode(0)->asEigenVector3d();

            // The start of the line element.
            auto const& a = bulk_element.getNode(0)->asEigenVector3d();

            // The end of the line element is the 2nd base node of the line,
            // which is the 2nd node.
            auto const& b = bulk_element.getNode(1)->asEigenVector3d();

            // x coincides either with the start or with the end of the line.
            // a + b - 2 * x evaluates to a - b or to b - a, respectively.
            // The formula assumes that the line is perfectly straight, even in
            // the LINE3 case.
            surface_element_normal = a + b - 2 * x;
        }
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
