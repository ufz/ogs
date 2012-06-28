/**
 * \file MathTools.cpp
 *
 *  Created on 2010-01-13 by Thomas Fischer
 */

#include "MathTools.h"

namespace MathLib {

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
	lambda = scpr<double,3> (u, v) / scpr<double,3> (v, v);

	// compute projected point
	double proj_pnt[3];
	for (size_t k(0); k<3; k++) proj_pnt[k] = a[k] + lambda * v[k];

	d0 = sqrt (sqrDist (proj_pnt, a));

	return sqrt (sqrDist (p, proj_pnt));
}

double sqrNrm2 (const GeoLib::Point* p0)
{
	return scpr<double,3> (p0->getCoords(), p0->getCoords());
}

double sqrDist (const GeoLib::Point* p0, const GeoLib::Point* p1)
{
	const double v[3] = {(*p1)[0] - (*p0)[0], (*p1)[1] - (*p0)[1], (*p1)[2] - (*p0)[2]};
	return scpr<double,3>(v,v);
}

double sqrDist(const double* p0, const double* p1)
{
	const double v[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
	return scpr<double,3>(v,v);
}

bool checkDistance(GeoLib::Point const &p0, GeoLib::Point const &p1, double squaredDistance)
{
	return (sqrDist(&p0, &p1) < squaredDistance);
}

float normalize(float min, float max, float val)
{
	return ((val-min)/static_cast<float>(max-min));
}

double getAngle (const double p0[3], const double p1[3], const double p2[3])
{
	const double v0[3] = {p0[0]-p1[0], p0[1]-p1[1], p0[2]-p1[2]};
	const double v1[3] = {p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]};

	// apply Cauchy Schwarz inequality
	return acos (scpr<double,3> (v0,v1) / (sqrt(scpr<double,3>(v0,v0)) * sqrt(scpr<double,3>(v1,v1))));
}

double calcTriangleArea(const double* p0, const double* p1, const double* p2)
{
	const double u0 (p2[0] - p0[0]);
	const double u1 (p2[1] - p0[1]);
	const double u2 (p2[2] - p0[2]);

	const double v0 (p1[0] - p0[0]);
	const double v1 (p1[1] - p0[1]);
	const double v2 (p1[2] - p0[2]);

	const double z0 (u1*v2 - u2*v1);
	const double z1 (u2*v0 - u0*v2);
	const double z2 (u0*v1 - u1*v0);

	return 0.5 * sqrt(z0*z0 + z1*z1 + z2 * z2);
}

double calcTetrahedronVolume(const double* x1, const double* x2, const double* x3, const double* x4)
{
	return fabs((x1[0] - x4[0]) * ((x2[1] - x4[1]) * (x3[2] - x4[2]) - (x2[2] - x4[2]) * (x3[1] - x4[1]))
	          - (x1[1] - x4[1]) * ((x2[0] - x4[0]) * (x3[2] - x4[2]) - (x2[2] - x4[2]) * (x3[0] - x4[0]))
	          + (x1[2] - x4[2]) * ((x2[0] - x4[0]) * (x3[1] - x4[1]) - (x2[1] - x4[1]) * (x3[0] - x4[0]))) / 6.0;
}

} // namespace
