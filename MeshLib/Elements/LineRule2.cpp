/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LineRule2.h"

#include "MathLib/MathTools.h"
#include "MeshLib/Node.h"

namespace MeshLib
{
const unsigned LineRule2::edge_nodes[1][2] = {
    {0, 1}  // Edge 0
};

double LineRule2::computeVolume(Node const* const* _nodes)
{
    return std::sqrt(MathLib::sqrDist(*_nodes[0], *_nodes[1]));
}

bool LineRule2::isPntInElement(Node const* const* nodes,
                               MathLib::Point3d const& pnt, double eps)
{
    auto const& a = *nodes[0];
    auto const& b = *nodes[1];

    if (MathLib::sqrDist(a, pnt) < eps * eps)
    {
        return true;
    }
    if (MathLib::sqrDist(b, pnt) < eps * eps)
    {
        return true;
    }

    double lambda;
    double distance_of_proj_pnt_to_a;
    auto const distance_from_line = MathLib::calcProjPntToLineAndDists(
        pnt, a, b, lambda, distance_of_proj_pnt_to_a);

    if (lambda >= 0 && lambda <= 1)
    {
        // pnt has been projected to the interior of the line
        return distance_from_line < eps;
    }
    else
    {
        // pnt has been projected outside the line segment
        // corner cases have already been treated in the beginning
        return false;
    }
}

unsigned LineRule2::identifyFace(Node const* const* _nodes,
                                 Node const* nodes[1])
{
    if (nodes[0] == _nodes[0])
    {
        return 0;
    }
    if (nodes[0] == _nodes[1])
    {
        return 1;
    }
    return std::numeric_limits<unsigned>::max();
}

ElementErrorCode LineRule2::validate(const Element* e)
{
    ElementErrorCode error_code;
    error_code[ElementErrorFlag::ZeroVolume] = hasZeroVolume(*e);
    return error_code;
}

}  // end namespace MeshLib
