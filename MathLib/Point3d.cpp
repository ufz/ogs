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

MathLib::Point3d operator*(Eigen::Matrix3d const& mat,
                           MathLib::Point3d const& p)
{
    auto const& result = (mat * p.asEigenVector3d()).eval();
    return MathLib::Point3d{{result[0], result[1], result[2]}};
}

double sqrDist(MathLib::Point3d const& p0, MathLib::Point3d const& p1)
{
    return (p0.asEigenVector3d() - p1.asEigenVector3d()).squaredNorm();
}

bool lessEq(Point3d const& a, Point3d const& b, double eps)
{
    auto absAndRelDiffLargerThanEps = [eps](double const u,
                                            double const v) -> bool
    {
        return std::abs(u - v) > eps * std::min(std::abs(v), std::abs(u)) &&
               std::abs(u - v) > eps;
    };

    return std::lexicographical_compare(
        a.data(), a.data() + 3, b.data(), b.data() + 3,
        [&absAndRelDiffLargerThanEps](auto const u, auto const v)
        {
            if (absAndRelDiffLargerThanEps(u, v))
            {
                return u <= v;
            }
            return true;
        });
}

extern const Point3d ORIGIN{{{0.0, 0.0, 0.0}}};
}  // namespace MathLib
