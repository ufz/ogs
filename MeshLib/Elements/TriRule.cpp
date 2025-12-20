// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "TriRule.h"

#include "MathLib/GeometricBasics.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
double TriRule::computeVolume(Node const* const* element_nodes)
{
    return MathLib::calcTriangleArea(*element_nodes[0], *element_nodes[1],
                                     *element_nodes[2]);
}

bool TriRule::isPntInElement(Node const* const* nodes,
                             MathLib::Point3d const& pnt, double eps)
{
    return MathLib::isPointInTriangle(pnt, *nodes[0], *nodes[1], *nodes[2],
                                      eps);
}

ElementErrorCode TriRule::validate(const Element* e)
{
    ElementErrorCode error_code;
    error_code[ElementErrorFlag::ZeroVolume] = hasZeroVolume(*e);
    error_code[ElementErrorFlag::NodeOrder] = !e->testElementNodeOrder();
    return error_code;
}

}  // end namespace MeshLib
