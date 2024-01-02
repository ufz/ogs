/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PrismRule.h"

#include "MathLib/GeometricBasics.h"
#include "MeshLib/Node.h"
#include "Quad.h"

namespace MeshLib
{
double PrismRule::computeVolume(Node const* const* element_nodes)
{
    return MathLib::calcTetrahedronVolume(*element_nodes[0],
                                          *element_nodes[1],
                                          *element_nodes[2],
                                          *element_nodes[3]) +
           MathLib::calcTetrahedronVolume(*element_nodes[1],
                                          *element_nodes[4],
                                          *element_nodes[2],
                                          *element_nodes[3]) +
           MathLib::calcTetrahedronVolume(*element_nodes[2],
                                          *element_nodes[4],
                                          *element_nodes[5],
                                          *element_nodes[3]);
}

bool PrismRule::isPntInElement(Node const* const* nodes,
                               MathLib::Point3d const& pnt,
                               double eps)
{
    return (MathLib::isPointInTetrahedron(
                pnt, *nodes[0], *nodes[1], *nodes[2], *nodes[3], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[1], *nodes[4], *nodes[2], *nodes[3], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[2], *nodes[4], *nodes[5], *nodes[3], eps));
}

ElementErrorCode PrismRule::validate(const Element* e)
{
    ElementErrorCode error_code;
    error_code[ElementErrorFlag::ZeroVolume] = hasZeroVolume(*e);

    for (unsigned i = 1; i < 4; ++i)
    {
        const auto* quad(dynamic_cast<const MeshLib::Quad*>(e->getFace(i)));
        if (quad)
        {
            error_code |= quad->validate();
        }
        else
        {
            error_code.set(ElementErrorFlag::NodeOrder);
        }
        delete quad;
    }
    error_code[ElementErrorFlag::NodeOrder] = !e->testElementNodeOrder();
    return error_code;
}

}  // end namespace MeshLib
