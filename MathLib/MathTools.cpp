/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MathTools.h"

#include <Eigen/Dense>
#include <cmath>

#include "Point3d.h"

namespace MathLib
{
double calcProjPntToLineAndDists(Point3d const& pp, Point3d const& pa,
                                 Point3d const& pb, double& lambda, double& d0)
{
    auto const& a = pa.asEigenVector3d();
    auto const& b = pb.asEigenVector3d();
    auto const& p = pp.asEigenVector3d();

    // g(lambda) = a + lambda v, v = b-a
    Eigen::Vector3d const v = b - a;

    // orthogonal projection: (p - g(lambda))^T * v = 0
    // <=> (a-p - lambda (b-a))^T * (b-a) = 0
    // <=> (a-p)^T * (b-a) = lambda (b-a)^T ) (b-a)
    lambda = (((p - a).transpose() * v) / v.squaredNorm())(0, 0);

    // compute projected point
    Eigen::Vector3d const proj_pnt = a + lambda * v;

    d0 = (proj_pnt - a).norm();

    return (p - proj_pnt).norm();
}

double getAngle(Point3d const& p0, Point3d const& p1, Point3d const& p2)
{
    auto const& b = p1.asEigenVector3d();
    Eigen::Vector3d const v0 = p0.asEigenVector3d() - b;
    Eigen::Vector3d const v1 = p2.asEigenVector3d() - b;

    // apply Cauchy Schwarz inequality
    return std::acos(v0.dot(v1) / (v0.norm() * v1.norm()));
}

}  // namespace MathLib
