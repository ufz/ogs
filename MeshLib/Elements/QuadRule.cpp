/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "QuadRule.h"

#include "MathLib/GeometricBasics.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
double QuadRule::computeVolume(Node const* const* element_nodes)
{
    return MathLib::calcTriangleArea(
               *element_nodes[0], *element_nodes[1], *element_nodes[2]) +
           MathLib::calcTriangleArea(
               *element_nodes[2], *element_nodes[3], *element_nodes[0]);
}

bool QuadRule::isPntInElement(Node const* const* nodes,
                              MathLib::Point3d const& pnt,
                              double eps)
{
    return (
        MathLib::isPointInTriangle(pnt, *nodes[0], *nodes[1], *nodes[2], eps) ||
        MathLib::isPointInTriangle(pnt, *nodes[0], *nodes[2], *nodes[3], eps));
}

ElementErrorCode QuadRule::validate(const Element* e)
{
    ElementErrorCode error_code;
    error_code[ElementErrorFlag::ZeroVolume] = hasZeroVolume(*e);
    Node const* const* element_nodes = e->getNodes();
    error_code[ElementErrorFlag::NonCoplanar] =
        (!MathLib::isCoplanar(*element_nodes[0],
                              *element_nodes[1],
                              *element_nodes[2],
                              *element_nodes[3]));
    // for collapsed quads (i.e. reduced to a line) this test might result
    // "false" as all four points are actually located on a line.
    if (!error_code[ElementErrorFlag::ZeroVolume])
    {
        error_code[ElementErrorFlag::NonConvex] =
            (!(MathLib::dividedByPlane(*element_nodes[0],
                                       *element_nodes[2],
                                       *element_nodes[1],
                                       *element_nodes[3]) &&
               MathLib::dividedByPlane(*element_nodes[1],
                                       *element_nodes[3],
                                       *element_nodes[0],
                                       *element_nodes[2])));
    }
    error_code[ElementErrorFlag::NodeOrder] = !e->testElementNodeOrder();
    return error_code;
}

}  // end namespace MeshLib
