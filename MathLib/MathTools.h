/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <Eigen/Eigen>
#include <cstddef>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace MathLib
{
template <typename T, std::size_t DIM> class TemplatePoint;
using Point3d = MathLib::TemplatePoint<double, 3>;

/**
 * standard inner product in R^N
 * \param v0 array of type T representing the vector
 * \param v1 array of type T representing the vector
 * */
template<typename T, int N> inline
T scalarProduct(T const * const v0, T const * const v1)
{
    T res (v0[0] * v1[0]);

#pragma omp parallel for reduction (+:res)
    for (int k = 1; k < N; k++)
    {
        res += v0[k] * v1[k];
    }
    return res;
}

template <> inline
double scalarProduct<double,3>(double const * const v0, double const * const v1)
{
    double res (v0[0] * v1[0]);
    for (std::size_t k(1); k < 3; k++)
    {
        res += v0[k] * v1[k];
    }
    return res;
}

template <typename T>
inline T scalarProduct(T const* const v0, T const* const v1, int const n)
{
    T res (v0[0] * v1[0]);

#pragma omp parallel for reduction (+:res)
    for (int k = 1; k < n; k++)
    {
        res += v0[k] * v1[k];
    }
    return res;
}

/**
 * calcProjPntToLineAndDists computes the orthogonal projection
 * of a point p to the line described by the points a and b,
 * \f$g(\lambda) = a + \lambda (b - a)\f$,
 * the distance between p and the projected point
 * and the distances between the projected point and the end
 * points pa, pb of the line
 * \param pp the (mesh) point
 * \param pa first point of line
 * \param pb second point of line
 * \param lambda the projected point described by the line equation above
 * \param d0 distance to the line point a
 * \returns the distance between pp and the orthogonal projection of pp
 */
double calcProjPntToLineAndDists(MathLib::Point3d const& pp,
                                 MathLib::Point3d const& pa,
                                 MathLib::Point3d const& pb, double& lambda,
                                 double& d0);

/**
 * Let \f$p_0, p_1, p_2 \in R^3\f$. The function getAngle
 * computes the angle between the edges \f$(p_0,p_1)\f$ and \f$(p_1,p_2)\f$
 * @param p0 start point of edge 0
 * @param p1 end point of edge 0 and start point of edge 1
 * @param p2 end point of edge 1
 * @return the angle between the edges
 */
double getAngle(Point3d const& p0, Point3d const& p1, Point3d const& p2);

/// Calculates the scalar triple (u x v) . w
double scalarTriple(Eigen::Vector3d const& u, Eigen::Vector3d const& v,
                    Eigen::Vector3d const& w);
}  // namespace MathLib
