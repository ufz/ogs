/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PyramidRule.h"

#include "MathLib/GeometricBasics.h"
#include "MeshLib/Node.h"
#include "Quad.h"

namespace MeshLib
{
double PyramidRule::computeVolume(Node const* const* element_nodes)
{
    return MathLib::calcTetrahedronVolume(*element_nodes[0],
                                          *element_nodes[1],
                                          *element_nodes[2],
                                          *element_nodes[4]) +
           MathLib::calcTetrahedronVolume(*element_nodes[2],
                                          *element_nodes[3],
                                          *element_nodes[0],
                                          *element_nodes[4]);
}

bool PyramidRule::isPntInElement(Node const* const* nodes,
                                 MathLib::Point3d const& pnt,
                                 double eps)
{
    return (MathLib::isPointInTetrahedron(
                pnt, *nodes[0], *nodes[1], *nodes[2], *nodes[4], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[0], *nodes[2], *nodes[3], *nodes[4], eps));
}

ElementErrorCode PyramidRule::validate(const Element* e)
{
    ElementErrorCode error_code;
    error_code[ElementErrorFlag::ZeroVolume] = hasZeroVolume(*e);

    std::unique_ptr<MeshLib::Quad const> const base{
        dynamic_cast<MeshLib::Quad const*>(e->getFace(4))};
    if (base)
    {
        error_code |= base->validate();
        error_code[ElementErrorFlag::NodeOrder] = !e->testElementNodeOrder();
    }
    else
    {
        error_code.set(ElementErrorFlag::NodeOrder);
    }

    return error_code;
}

}  // end namespace MeshLib
