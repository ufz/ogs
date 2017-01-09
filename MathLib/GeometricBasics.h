/**
 *
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
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

enum TriangleTest
{
    GAUSS, BARYCENTRIC
};

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

/**
 * Tests if the given point p is within the triangle, defined by its edge nodes
 * a, b and c.Using two eps-values it is possible to test an 'epsilon'
 * neighbourhood around the triangle as well as an 'epsilon' outside the
 * triangles plane.
 * @param p test point
 * @param a edge node of triangle
 * @param b edge node of triangle
 * @param c edge node of triangle
 * @param eps_pnt_out_of_plane eps allowing for p to be slightly off the plane
 * spanned by abc
 * @param eps_pnt_out_of_tri eps allowing for p to be slightly off outside of
 * abc
 * @param algorithm defines the method to use
 * @return true if the test point p is within the 'epsilon'-neighbourhood of the
 * triangle
 */
bool isPointInTriangle(
    MathLib::Point3d const& p,
    MathLib::Point3d const& a,
    MathLib::Point3d const& b,
    MathLib::Point3d const& c,
    double eps_pnt_out_of_plane = std::numeric_limits<float>::epsilon(),
    double eps_pnt_out_of_tri = std::numeric_limits<float>::epsilon(),
    MathLib::TriangleTest algorithm = MathLib::GAUSS);

/**
 * Tests if the given point p is within the triangle, defined by its edge nodes
 * a, b and c.
 * Using two eps-values it is possible to test an 'epsilon' neighbourhood around
 * the triangle
 * as well as an 'epsilon' outside the triangles plane.
 * @param p test point
 * @param a edge node of triangle
 * @param b edge node of triangle
 * @param c edge node of triangle
 * @param eps_pnt_out_of_plane eps allowing for p to be slightly off the plane
 * spanned by abc ((orthogonal distance to the plane spaned by triangle)
 * @param eps_pnt_out_of_tri eps allowing for p to be slightly off outside of
 * abc
 * @return true if the test point p is within the 'epsilon'-neighbourhood of the
 * triangle
 */
bool gaussPointInTriangle(
    MathLib::Point3d const& p, MathLib::Point3d const& a,
    MathLib::Point3d const& b, MathLib::Point3d const& c,
    double eps_pnt_out_of_plane = std::numeric_limits<float>::epsilon(),
    double eps_pnt_out_of_tri = std::numeric_limits<float>::epsilon());

/**
 * Tests if the given point p is within the triangle, defined by its edge nodes
 * a, b and c.
 * Using two eps-values it is possible to test an 'epsilon' neighbourhood around
 * the triangle
 * as well as an 'epsilon' outside the triangles plane.
 * Algorithm based on "Fundamentals of Computer Graphics" by Peter Shirley.
 * @param p test point
 * @param a edge node of triangle
 * @param b edge node of triangle
 * @param c edge node of triangle
 * @param eps_pnt_out_of_plane eps allowing for p to be slightly off the plane
 * spanned by abc
 * @param eps_pnt_out_of_tri eps allowing for p to be slightly off outside of
 * abc
 * @return true if the test point p is within the 'epsilon'-neighbourhood of the
 * triangle
 */
bool barycentricPointInTriangle(
    MathLib::Point3d const& p, MathLib::Point3d const& a,
    MathLib::Point3d const& b, MathLib::Point3d const& c,
    double eps_pnt_out_of_plane = std::numeric_limits<float>::epsilon(),
    double eps_pnt_out_of_tri = std::numeric_limits<float>::epsilon());

/// Checks if the point \f$p'\f$ is in the triangle defined by the points
/// \f$a', b', c'\f$, where the \f$p', a', b', c' \f$ are the orthogonal
/// projections to the \f$x\f$-\f$y\f$ plane of the points \f$p, a, b, c\f$,
/// respectively.
bool isPointInTriangleXY(MathLib::Point3d const& p, MathLib::Point3d const& a,
MathLib::Point3d const& b, MathLib::Point3d const& c);

/**
 * Checks if a and b can be placed on a plane such that c and d lie on different
 * sides of said plane. (In 2D space this checks if c and d are on different
 * sides of a line through a and b.)
 * @param a first point on plane
 * @param b second point on plane
 * @param c first point to test
 * @param d second point to test
 * @return true, if such a plane can be found, false otherwise
 */
bool dividedByPlane(const MathLib::Point3d& a, const MathLib::Point3d& b,
                    const MathLib::Point3d& c, const MathLib::Point3d& d);

/// Checks if the four given points are located on a plane.
bool isCoplanar(const MathLib::Point3d& a, const MathLib::Point3d& b,
                const MathLib::Point3d& c, const MathLib::Point3d& d);

}  // end namespace MathLib

#endif
