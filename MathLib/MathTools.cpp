/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <cmath>

#include "MathTools.h"

namespace MathLib
{

double calcProjPntToLineAndDists(const double p[3], const double a[3],
        const double b[3], double &lambda, double &d0)
{
    // g (lambda) = a + lambda v, v = b-a
    double v[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    // orthogonal projection: (g(lambda)-p) * v = 0 => in order to compute lambda we define a help vector u
    double u[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};
    lambda = scalarProduct<double,3> (u, v) / scalarProduct<double,3> (v, v);

    // compute projected point
    double proj_pnt[3];
    for (std::size_t k(0); k < 3; k++)
        proj_pnt[k] = a[k] + lambda * v[k];

    d0 = std::sqrt (sqrDist (proj_pnt, a));

    return std::sqrt (sqrDist (p, proj_pnt));
}

double getAngle (const double p0[3], const double p1[3], const double p2[3])
{
    const double v0[3] = {p0[0]-p1[0], p0[1]-p1[1], p0[2]-p1[2]};
    const double v1[3] = {p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]};

    // apply Cauchy Schwarz inequality
    return std::acos (scalarProduct<double,3> (v0,v1) / (std::sqrt(scalarProduct<double,3>(v0,v0)) * sqrt(scalarProduct<double,3>(v1,v1))));
}



} // namespace
