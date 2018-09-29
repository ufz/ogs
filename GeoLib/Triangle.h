/**
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

#include "MathLib/Point3d.h"

namespace GeoLib {

class Point;

/** \brief Class Triangle consists of a reference to a point vector and
 * a vector that stores the indices in the point vector.
 * A surface is composed by triangles. The class Surface stores the position
 * of pointers to the points of triangles in the _pnt_ids vector.
 * */
class Triangle final
{
public:
    /**
     * construction of object, initialization of reference to point vector,
     * saves the three indices describing a triangle
     */
    Triangle(std::vector<Point*> const& pnt_vec,
             std::size_t pnt_a,
             std::size_t pnt_b,
             std::size_t pnt_c);

    /**
     * saves three indices describing a triangle
     * */
    void setTriangle (std::size_t pnt_a, std::size_t pnt_b, std::size_t pnt_c);

    /** \brief const access operator to access the index
     * of the i-th triangle point
    */
    const std::size_t& operator[](std::size_t i) const
    {
        assert (i < 3);
        return _pnt_ids[i];
    }

    /**
     * \brief const access operator to access the i-th triangle Point
     */
    const Point* getPoint(std::size_t i) const
    {
        assert (i < 3);
        return _pnts[_pnt_ids[i]];
    }

    /**
     * Checks if point q is within the triangle, uses GeoLib::isPointInTriangle().
     * @param q The input point.
     * @param eps Checks the 'epsilon'-neighbourhood
     * @return true, if point is in triangle, else false
     */
    bool containsPoint(
        MathLib::Point3d const& q,
        double eps = std::numeric_limits<float>::epsilon()) const;

    /**
     * projects the triangle points to the x-y-plane and
     * checks if point pnt is contained into the triangle
     * @param pnt the point to test for
     * @return true, if the point is into the projected triangle
     */
    bool containsPoint2D(Point const& pnt) const;

private:
    /// a vector of pointers to points the triangle is based on
    std::vector<Point*> const& _pnts;
    /// position of pointers to the geometric points
    std::array<std::size_t, 3> _pnt_ids;
};

void getPlaneCoefficients(Triangle const& tri, double c[3]);

} // end namespace GeoLib
