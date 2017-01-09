/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "TriRule3.h"

#include <logog/include/logog.hpp>

#include "MathLib/GeometricBasics.h"

#include "MeshLib/Node.h"

namespace MeshLib {

const unsigned TriRule3::edge_nodes[3][2] =
{
        {0, 1}, // Edge 0
        {1, 2}, // Edge 1
        {2, 0}, // Edge 2
};

double TriRule3::computeVolume(Node const* const* _nodes)
{
    return MathLib::calcTriangleArea(*_nodes[0], *_nodes[1], *_nodes[2]);
}

bool TriRule3::isPntInElement(Node const* const* nodes,
                              MathLib::Point3d const& pnt, double eps)
{
    return MathLib::isPointInTriangle(pnt, *nodes[0], *nodes[1], *nodes[2],
                                      eps);
}

unsigned TriRule3::identifyFace(Node const* const* _nodes, Node* nodes[3])
{
    for (unsigned i=0; i<3; i++)
    {
        unsigned flag(0);
        for (unsigned j=0; j<2; j++)
            for (unsigned k=0; k<2; k++)
                if (_nodes[edge_nodes[i][j]] == nodes[k])
                    flag++;
        if (flag==2)
            return i;
    }
    return std::numeric_limits<unsigned>::max();
}

ElementErrorCode TriRule3::validate(const Element* e)
{
    ElementErrorCode error_code;
    error_code[ElementErrorFlag::ZeroVolume] = e->hasZeroVolume();
    error_code[ElementErrorFlag::NodeOrder]  = !e->testElementNodeOrder();
    return error_code;
}

} // end namespace MeshLib
