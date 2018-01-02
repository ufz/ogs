/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PointRule1.h"

#include "MathLib/Point3d.h"
#include "MeshLib/Node.h"

namespace MeshLib {

const unsigned PointRule1::edge_nodes[1][1] =
{
    {0}
};

double PointRule1::computeVolume(Node const* const* /*_nodes*/)
{
    return 0;
}

bool PointRule1::isPntInElement(Node const* const* nodes,
                                MathLib::Point3d const& pnt, double eps)
{
    double const dist = MathLib::sqrDist(*nodes[0], pnt);
    return (dist < eps);
}

unsigned PointRule1::identifyFace(Node const* const* _nodes, Node* nodes[1])
{
    if (nodes[0] == _nodes[0])
        return 0;
    return std::numeric_limits<unsigned>::max();
}

ElementErrorCode PointRule1::validate(const Element* e)
{
    ElementErrorCode error_code;
    error_code[ElementErrorFlag::ZeroVolume] = e->hasZeroVolume();
    return error_code;
}

} // end namespace MeshLib
