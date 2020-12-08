/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    auto const a =
        Eigen::Map<Eigen::Vector3d>(const_cast<double*>(pa.getCoords()));
    auto const b =
        Eigen::Map<Eigen::Vector3d>(const_cast<double*>(pb.getCoords()));
    auto const p =
        Eigen::Map<Eigen::Vector3d>(const_cast<double*>(pp.getCoords()));
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

double getAngle (const double p0[3], const double p1[3], const double p2[3])
{
    const double v0[3] = {p0[0]-p1[0], p0[1]-p1[1], p0[2]-p1[2]};
    const double v1[3] = {p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]};

    // apply Cauchy Schwarz inequality
    return std::acos (scalarProduct<double,3> (v0,v1) / (std::sqrt(scalarProduct<double,3>(v0,v0)) * sqrt(scalarProduct<double,3>(v1,v1))));
}

double scalarTriple(Eigen::Vector3d const& u, Eigen::Vector3d const& v,
                    Eigen::Vector3d const& w)
{
    return u.cross(v).dot(w);
}

}  // namespace MathLib
