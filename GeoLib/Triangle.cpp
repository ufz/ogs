/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Triangle.h"

#include <Eigen/Dense>

#include "Point.h"
#include "AnalyticalGeometry.h"

#include "MathLib/GeometricBasics.h"

namespace GeoLib {

Triangle::Triangle (std::vector<Point *> const &pnt_vec,
    std::size_t pnt_a, std::size_t pnt_b, std::size_t pnt_c) :
    pnts_(pnt_vec), pnt_ids_( {{pnt_a, pnt_b, pnt_c}} )
{
    assert(!pnts_.empty());
    assert (pnt_a < pnts_.size() && pnt_b < pnts_.size() && pnt_c < pnts_.size());
}

bool Triangle::containsPoint(MathLib::Point3d const& q, double eps) const
{
    GeoLib::Point const& a(*(pnts_[pnt_ids_[0]]));
    GeoLib::Point const& b(*(pnts_[pnt_ids_[1]]));
    GeoLib::Point const& c(*(pnts_[pnt_ids_[2]]));
    return MathLib::isPointInTriangle(q, a, b, c, eps);
}
}  // namespace GeoLib
