/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TetRule.h"

#include "MathLib/GeometricBasics.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
double TetRule::computeVolume(Node const* const* element_nodes)
{
    return MathLib::calcTetrahedronVolume(*element_nodes[0], *element_nodes[1],
                                          *element_nodes[2], *element_nodes[3]);
}

bool TetRule::isPntInElement(Node const* const* nodes,
                             MathLib::Point3d const& pnt, double eps)
{
    return MathLib::isPointInTetrahedron(pnt, *nodes[0], *nodes[1], *nodes[2],
                                         *nodes[3], eps);
}

ElementErrorCode TetRule::validate(const Element* e)
{
    ElementErrorCode error_code;
    error_code[ElementErrorFlag::ZeroVolume] = hasZeroVolume(*e);
    error_code[ElementErrorFlag::NodeOrder] = !e->testElementNodeOrder();
    return error_code;
}

}  // end namespace MeshLib
