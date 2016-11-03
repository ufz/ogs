/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MATHTOOLS_H_
#define MATHTOOLS_H_

#include <cstddef>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef Q_MOC_RUN // to avoid Qt4 bug, https://bugreports.qt.io/browse/QTBUG-22829
#include <boost/math/constants/constants.hpp>
#endif

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

/** squared dist between double arrays p0 and p1 (size of arrays is 3) */
inline
double sqrDist(const double* p0, const double* p1)
{
    const double v[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
    return scalarProduct<double,3>(v,v);
}

/**
 * Let \f$p_0, p_1, p_2 \in R^3\f$. The function getAngle
 * computes the angle between the edges \f$(p_0,p_1)\f$ and \f$(p_1,p_2)\f$
 * @param p0 start point of edge 0
 * @param p1 end point of edge 0 and start point of edge 1
 * @param p2 end point of edge 1
 * @return the angle between the edges
 */
double getAngle (const double p0[3], const double p1[3], const double p2[3]);

/// converts the given degrees to radians
inline double to_radians(double degrees) {
    return degrees*boost::math::constants::pi<double>()/180.;
}

template<typename Type> Type limitValueInInterval(const Type variable,
                                                  const Type lower_bound,
                                                  const Type upper_bound)
{
    if (variable < lower_bound)
        return lower_bound;
    if (variable > upper_bound)
        return upper_bound;
    return variable;
}
} // namespace

#endif /* MATHTOOLS_H_ */
