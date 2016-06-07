/**
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef GEOMETRIC_BASICS_H_
#define GEOMETRIC_BASICS_H_

#include <cstddef>

namespace MathLib
{

template <typename T, std::size_t DIM> class TemplatePoint;
typedef MathLib::TemplatePoint<double,3> Point3d;

/**
 * Checks if a point p is on the left or right side of a plane spanned by three
 * points a, b, c.
 * @param p point to test
 * @param a first point on plane
 * @param b second point on plane
 * @param c third point on plane
 * @return If the triangle abc is ordered counterclockwise when viewed from p,
 * the method will return a negative value,
 * otherwise it will return a positive value. If p is coplanar with abc, the
 * function will return 0.
 */
double orientation3d(MathLib::Point3d const& p,
                     MathLib::Point3d const& a,
                     MathLib::Point3d const& b,
                     MathLib::Point3d const& c);

/**
 * Calculates the volume of a tetrahedron.
 * The formula is V=1/6*|a(b x c)| with a=x1->x2, b=x1->x3 and c=x1->x4.
 */
double calcTetrahedronVolume(MathLib::Point3d const& x1,
                             MathLib::Point3d const& x2,
                             MathLib::Point3d const& x3,
                             MathLib::Point3d const& x4);

}  // end namespace MathLib

#endif
