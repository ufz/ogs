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
#include <limits>

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

/**
 * Calculates the area of the triangle defined by its edge nodes a, b and c.
 * The formula is \f$A= \frac{1}{2} \cdot |u \times v|\f$, i.e. half of the area of the
 * parallelogram specified by the vectors\f$u=b-a\f$ and \f$v=c-a\f$.
 */
double calcTriangleArea(MathLib::Point3d const& a, MathLib::Point3d const& b,
                        MathLib::Point3d const& c);

/**
 * Tests if the given point p is located within a tetrahedron spanned by points
 * a, b, c, d.
 * If the tet specified by a, b, c, d is degenerated (i.e. all points are
 * coplanar) the function
 * will return false because there is no meaningful concept of "inside" for such
 * elements.
 * @param p test point
 * @param a edge node of tetrahedron
 * @param b edge node of tetrahedron
 * @param c edge node of tetrahedron
 * @param d edge node of tetrahedron
 * @param eps Accounts for numerical inaccuracies by allowing a point to be
 * slightly outside of the element and still be regarded as inside.
 * @return true if the test point p is not located outside of abcd (i.e. inside
 * or on a plane/edge).
 */
bool isPointInTetrahedron(MathLib::Point3d const& p, MathLib::Point3d const& a,
                          MathLib::Point3d const& b, MathLib::Point3d const& c,
                          MathLib::Point3d const& d,
                          double eps = std::numeric_limits<double>::epsilon());

}  // end namespace MathLib

#endif
