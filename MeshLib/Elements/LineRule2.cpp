/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LineRule2.h"

#include <logog/include/logog.hpp>

#include "MeshLib/Node.h"

namespace MeshLib {

const unsigned LineRule2::edge_nodes[1][2] =
{
    {0, 1} // Edge 0
};

double LineRule2::computeVolume(Node const* const* _nodes)
{
    return sqrt(MathLib::sqrDist(_nodes[0]->getCoords(), _nodes[1]->getCoords()));
}

bool LineRule2::isPntInElement(Node const* const* nodes,
                               MathLib::Point3d const& pnt, double eps)
{
    double tmp;
    double tmp_dst(0);
    double const dist = MathLib::calcProjPntToLineAndDists(
        pnt.getCoords(), nodes[0]->getCoords(), nodes[1]->getCoords(), tmp,
        tmp_dst);
    return (dist < eps);
}

unsigned LineRule2::identifyFace(Node const* const* _nodes, Node* nodes[1])
{
    if (nodes[0] == _nodes[0])
        return 0;
    if (nodes[0] == _nodes[1])
        return 1;
    return std::numeric_limits<unsigned>::max();
}

ElementErrorCode LineRule2::validate(const Element* e)
{
    ElementErrorCode error_code;
    error_code[ElementErrorFlag::ZeroVolume] = e->hasZeroVolume();
    return error_code;
}

} // end namespace MeshLib
