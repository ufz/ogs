/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#include "Tests/GeoLib/AutoCheckGenerators.h"

namespace autocheck
{

// reflect point p on the point c in x-y plane
MathLib::Point3d reflect(MathLib::Point3d const& c, MathLib::Point3d const& p)
{
    return MathLib::Point3d(
        std::array<double, 3>{{2 * c[0] - p[0], 2 * c[1] - p[1], 0.0}});
}

GeoLib::LineSegment translate(MathLib::Vector3 const& translation,
                              GeoLib::LineSegment const& line_seg)
{
    auto a = std::make_unique<GeoLib::Point>(line_seg.getBeginPoint());
    auto b = std::make_unique<GeoLib::Point>(line_seg.getEndPoint());
    for (std::size_t k(0); k<3; ++k)
        (*a)[k] += translation[k];
    for (std::size_t k(0); k<3; ++k)
        (*b)[k] += translation[k];
    return GeoLib::LineSegment{a.release(), b.release(), true};
}

}  // namespace autocheck
