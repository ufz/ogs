/**
 * \file
 * \author Thomas Fischer
 * \date   2010-01-13
 * \brief  Implementation of math helper functions.
 *
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MathTools.h"

namespace MathLib
{
void crossProd(const double u[3], const double v[3], double r[3])
{
	r[0] = u[1] * v[2] - u[2] * v[1];
	r[1] = u[2] * v[0] - u[0] * v[2];
	r[2] = u[0] * v[1] - u[1] * v[0];
}

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
	for (size_t k(0); k < 3; k++)
		proj_pnt[k] = a[k] + lambda * v[k];

	d0 = sqrt (sqrDist (proj_pnt, a));

	return sqrt (sqrDist (p, proj_pnt));
}

double sqrDist(const double* p0, const double* p1)
{
	const double v[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
	return scalarProduct<double,3>(v,v);
}

float normalize(float min, float max, float val)
{
	return (val - min) / static_cast<float>(max - min);
}

double getAngle (const double p0[3], const double p1[3], const double p2[3])
{
	const double v0[3] = {p0[0]-p1[0], p0[1]-p1[1], p0[2]-p1[2]};
	const double v1[3] = {p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]};

	// apply Cauchy Schwarz inequality
	return acos (scalarProduct<double,3> (v0,v1) / (sqrt(scalarProduct<double,3>(v0,v0)) * sqrt(scalarProduct<double,3>(v1,v1))));
}

} // namespace
