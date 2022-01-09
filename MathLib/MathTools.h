/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <cstddef>

namespace MathLib
{
template <typename T, std::size_t DIM> class TemplatePoint;
using Point3d = MathLib::TemplatePoint<double, 3>;

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
}  // namespace MathLib
