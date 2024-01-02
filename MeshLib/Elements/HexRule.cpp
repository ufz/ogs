/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "HexRule.h"

#include "MathLib/GeometricBasics.h"
#include "MeshLib/Node.h"
#include "Quad.h"

namespace MeshLib
{
double HexRule::computeVolume(Node const* const* element_nodes)
{
    return MathLib::calcTetrahedronVolume(*element_nodes[4],
                                          *element_nodes[7],
                                          *element_nodes[5],
                                          *element_nodes[0]) +
           MathLib::calcTetrahedronVolume(*element_nodes[5],
                                          *element_nodes[3],
                                          *element_nodes[1],
                                          *element_nodes[0]) +
           MathLib::calcTetrahedronVolume(*element_nodes[5],
                                          *element_nodes[7],
                                          *element_nodes[3],
                                          *element_nodes[0]) +
           MathLib::calcTetrahedronVolume(*element_nodes[5],
                                          *element_nodes[7],
                                          *element_nodes[6],
                                          *element_nodes[2]) +
           MathLib::calcTetrahedronVolume(*element_nodes[1],
                                          *element_nodes[3],
                                          *element_nodes[5],
                                          *element_nodes[2]) +
           MathLib::calcTetrahedronVolume(*element_nodes[3],
                                          *element_nodes[7],
                                          *element_nodes[5],
                                          *element_nodes[2]);
}

bool HexRule::isPntInElement(Node const* const* nodes,
                             MathLib::Point3d const& pnt,
                             double eps)
{
    return (MathLib::isPointInTetrahedron(
                pnt, *nodes[4], *nodes[7], *nodes[5], *nodes[0], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[5], *nodes[3], *nodes[1], *nodes[0], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[5], *nodes[7], *nodes[3], *nodes[0], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[5], *nodes[7], *nodes[6], *nodes[2], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[1], *nodes[3], *nodes[5], *nodes[2], eps) ||
            MathLib::isPointInTetrahedron(
                pnt, *nodes[3], *nodes[7], *nodes[5], *nodes[2], eps));
}

ElementErrorCode HexRule::validate(const Element* e)
{
    ElementErrorCode error_code;
    error_code[ElementErrorFlag::ZeroVolume] = hasZeroVolume(*e);

    for (unsigned i = 0; i < 6; ++i)
    {
        if (error_code.all())
        {
            break;
        }

        const MeshLib::Element* quad(e->getFace(i));
        error_code |= quad->validate();
        delete quad;
    }
    error_code[ElementErrorFlag::NodeOrder] = !e->testElementNodeOrder();
    return error_code;
}

}  // end namespace MeshLib
