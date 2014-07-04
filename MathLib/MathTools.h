/**
 * \file
 * \author Thomas Fischer
 * \date   2010-01-13
 * \brief  Definition of math helper functions.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHTOOLS_H_
#define MATHTOOLS_H_

#include <cmath>
#include <limits>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "TemplatePoint.h"

namespace MathLib
{
/**
 * standard inner product in R^N
 * \param v0 array of type T representing the vector
 * \param v1 array of type T representing the vector
 * */
template<typename T, int N> inline
T scalarProduct(T const * const v0, T const * const v1)
{
	T res (v0[0] * v1[0]);
#ifdef _OPENMP
	OPENMP_LOOP_TYPE k;

#pragma omp parallel for reduction (+:res)
	for (k = 1; k<N; k++) {
		res += v0[k] * v1[k];
	}
#else
	for (std::size_t k(1); k < N; k++)
		res += v0[k] * v1[k];
#endif
	return res;
}

template <> inline
double scalarProduct<double,3>(double const * const v0, double const * const v1)
{
	double res (v0[0] * v1[0]);
	for (std::size_t k(1); k < 3; k++)
		res += v0[k] * v1[k];
	return res;
}

template<typename T> inline
T scalarProduct(T const * const v0, T const * const v1, unsigned n)
{
	T res (v0[0] * v1[0]);
#ifdef _OPENMP
	OPENMP_LOOP_TYPE k;

#pragma omp parallel for reduction (+:res)
#ifdef WIN32
#pragma warning ( push )
#pragma warning ( disable: 4018 )
#endif
	for (k = 1; k<n; k++) {
		res += v0[k] * v1[k];
	}
#ifdef WIN32
#pragma warning ( pop )
#endif
#else
	for (std::size_t k(1); k < n; k++)
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

template <typename POINT_T>
typename POINT_T::FP_T sqrDist(POINT_T const& p0, POINT_T const& p1)
{
	typename POINT_T::FP_T const v[3] = {p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]};
	return MathLib::scalarProduct<typename POINT_T::FP_T,3>(v,v);
}

template <typename T, std::size_t DIM>
bool operator==(TemplatePoint<T,DIM> const& a, TemplatePoint<T,DIM> const& b)
{
	T const sqr_dist(sqrDist(a,b));
	return (sqr_dist < pow(std::numeric_limits<T>::epsilon(),2));
}

/** squared dist between double arrays p0 and p1 (size of arrays is 3) */
double sqrDist(const double* p0, const double* p1);

/** Distance between points p0 and p1 in the maximum norm. */
template <typename T>
T maxNormDist(const MathLib::TemplatePoint<T>* p0, const MathLib::TemplatePoint<T>* p1)
{
	const T x = fabs((*p1)[0] - (*p0)[0]);
	const T y = fabs((*p1)[1] - (*p0)[1]);
	const T z = fabs((*p1)[2] - (*p0)[2]);

	return std::max(x, std::max(y, z));
}

/** linear normalisation of val from [min, max] into [0,1] */
float normalize(float min, float max, float val);

/**
 * Let \f$p_0, p_1, p_2 \in R^3\f$. The function getAngle
 * computes the angle between the edges \f$(p_0,p_1)\f$ and \f$(p_1,p_2)\f$
 * @param p0 start point of edge 0
 * @param p1 end point of edge 0 and start point of edge 1
 * @param p2 end point of edge 1
 * @return the angle between the edges
 */
double getAngle (const double p0[3], const double p1[3], const double p2[3]);


/**
 * simple power function that takes as a second argument an integer instead of a float
 * @param base basis of the expression
 * @param exp exponent of the expression
 * @return base^exp
 */
template <typename T> inline
T fastpow (T base, std::size_t exp)
{
	T result (base);
	if (exp == 0)
		result = static_cast<T>(1);
	for (std::size_t k(1); k < exp; k++)
		result *= base;
	return result;
}

} // namespace

#endif /* MATHTOOLS_H_ */
