/**
 * \file MathTools.h
 *
 *  Created on 2010-01-13 by Thomas Fischer
 */

#ifndef MATHTOOLS_H_
#define MATHTOOLS_H_

#include <vector>
#include <cmath>
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "Point.h"

namespace MathLib {

/**
 * standard inner product in R^N
 * \param v0 array of type T representing the vector
 * \param v1 array of type T representing the vector
 * */
template<typename T, int N> inline
T scpr(T const * const v0, T const * const v1)
{
	T res (v0[0]*v1[0]);
#ifdef _OPENMP
	OPENMP_LOOP_TYPE k;

	#pragma omp parallel for reduction (+:res)
	for (k = 1; k<N; k++) {
		res += v0[k] * v1[k];
	}
#else
	for (size_t k(1); k < N; k++)
		res += v0[k] * v1[k];
#endif
	return res;
}

template <> inline
double scpr<double,3>(double const * const v0, double const * const v1)
{
	double res (v0[0]*v1[0]);
	for (size_t k(1); k < 3; k++)
		res += v0[k] * v1[k];
	return res;
}

template<typename T> inline
T scpr(T const * const v0, T const * const v1, unsigned n)
{
	T res (v0[0]*v1[0]);
#ifdef _OPENMP
	OPENMP_LOOP_TYPE k;

	#pragma omp parallel for reduction (+:res)
	for (k = 1; k<n; k++) {
		res += v0[k] * v1[k];
	}
#else
	for (size_t k(1); k < n; k++)
		res += v0[k] * v1[k];
#endif
	return res;
}


/**
 * computes the cross (or vector) product of the 3d vectors u and v
 * the result is given in the vector r
 */
void crossProd (const double u[3], const double v[3], double r[3]);

/**
 * calcProjPntToLineAndDists computes the orthogonal projection
 * of a point p to the line described by the points a and b,
 * \f$g(\lambda) = a + \lambda (b - a)\f$,
 * the distance between p and the projected point
 * and the distances between the projected point and the end
 * points a, b of the line
 * \param p the (mesh) point
 * \param a first point of line
 * \param b second point of line
 * \param lambda the projected point described by the line equation above
 * \param d0 distance to the line point a
 * \returns the distance between p and the orthogonal projection of p
 */
double calcProjPntToLineAndDists(const double p[3], const double a[3],
		const double b[3], double &lambda, double &d0);

/**
 * Checks if two points are within a given distance of each other
 * @param p0 The first point
 * @param p1 the second point
 * @param squaredDistance The square of the distance within which the two points should be
 * @return true if p1 and p2 are within the given distance of each other, false otherwise
 */
bool checkDistance(GeoLib::Point const &p0, GeoLib::Point const &p1, double squaredDistance);

/** squared euklid norm of the vector p0 */
double sqrNrm2(const GeoLib::Point* const p0);

/** squared dist between GeoLib::Points p0 and p1 */
double sqrDist(const GeoLib::Point* p0, const GeoLib::Point* p1);

/** squared dist between double arrays p0 and p1 (size of arrays is 3) */
double sqrDist(const double* p0, const double* p1);

/** linear normalisation of val from [min, max] into [0,1] */
float normalize(float min, float max, float val);

/**
 * Let \f$p_0, p_1, p_2 \in \mathbb{R}^3\f$. The function getAngle
 * computes the angle between the edges \f$(p_0,p_1)\f$ and \f$(p_1,p_2)\f$
 * @param p0 start point of edge 0
 * @param p1 end point of edge 0 and start point of edge 1
 * @param p2 end point of edge 1
 * @return the angle between the edges
 */
double getAngle (const double p0[3], const double p1[3], const double p2[3]);

/**
 * Calculates the area of a triangle.
 * The formula is A=.5*|u x v|, i.e. half of the area of the parallelogram specified by u=p0->p1 and v=p0->p2.
 */
double calcTriangleArea(const double* p0, const double* p1,	const double* p2);

/**
 * Calculates the volume of a tetrahedron.
 * The formula is V=1/6*|a(b x c)| with a=x1->x2, b=x1->x3 and c=x1->x4.
 */
double calcTetrahedronVolume(const double* x1, const double* x2, const double* x3, const double* x4);

/**
 * simple power function that takes as a second argument an integer instead of a float
 * @param base basis of the expression
 * @param exp exponent of the expression
 * @return base^exp
 */
template <typename T> inline
T fastpow (T base, size_t exp)
{
	T result (base);
	if (exp == 0) result = static_cast<T>(1);
	for (size_t k(1); k<exp; k++) {
		result *= base;
	}
	return result;
}

} // namespace

#endif /* MATHTOOLS_H_ */
