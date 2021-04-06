/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Triangle.h"

#include <Eigen/Dense>

#include "AnalyticalGeometry.h"
#include "MathLib/GeometricBasics.h"
#include "Point.h"

namespace GeoLib
{
Triangle::Triangle(std::vector<Point*> const& pnt_vec, std::size_t pnt_a,
                   std::size_t pnt_b, std::size_t pnt_c)
    : _pnts(pnt_vec), _pnt_ids({{pnt_a, pnt_b, pnt_c}})
{
    assert(!_pnts.empty());
    assert(pnt_a < _pnts.size() && pnt_b < _pnts.size() &&
           pnt_c < _pnts.size());
}

bool Triangle::containsPoint(MathLib::Point3d const& q, double eps) const
{
    GeoLib::Point const& a(*(_pnts[_pnt_ids[0]]));
    GeoLib::Point const& b(*(_pnts[_pnt_ids[1]]));
    GeoLib::Point const& c(*(_pnts[_pnt_ids[2]]));
    return MathLib::isPointInTriangle(q, a, b, c, eps);
}
}  // namespace GeoLib
