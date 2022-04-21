/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Point3d.h"

namespace MathLib
{
Point3d::Point3d() : x_({0, 0, 0}) {}

Point3d::Point3d(std::array<double, 3> x) : x_(x[0], x[1], x[2]) {}

double sqrDist(MathLib::Point3d const& p0, MathLib::Point3d const& p1)
{
    return (p0.asEigenVector3d() - p1.asEigenVector3d()).squaredNorm();
}

extern const Point3d ORIGIN{{{0.0, 0.0, 0.0}}};
}  // namespace MathLib
