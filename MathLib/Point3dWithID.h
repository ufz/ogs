/**
 * \file
 * \date   2015-05-18
 * \brief  Definition of the Point3d class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <limits>

#include "Point3d.h"

namespace MathLib
{
/**
 * Class Point3dWithID is derived from class Point3d in
 * order to extend the class Point3d with an id.
 */
class Point3dWithID: public Point3d {
public:
    /// Constructs a point with the coordinates x0, x1 and x2 and the provided
    /// id.
    /// @param x0 x coordinate of point
    /// @param x1 y coordinate of point
    /// @param x2 z coordinate of point
    /// @param id the id of the object [default: max of std::size_t]
    Point3dWithID(double x0, double x1, double x2,
        std::size_t id = std::numeric_limits<std::size_t>::max())
        : Point3d(std::array<double,3>({{x0, x1, x2}})), _id(id)
    {}

    /// Constructs a point using std::array<double,3> as coordinates and
    /// the provided id.
    /// @param coords coordinates of the point
    /// @param id the id of the object [default: max of std::size_t]
    Point3dWithID(std::array<double,3> const& coords,
        std::size_t id = std::numeric_limits<std::size_t>::max())
        : Point3d(coords), _id(id)
    {}

    /// Constructs a point with the same coordinates as the given
    /// Point3d pnt and the provided id.
    /// @param pnt a MathLib::Point3d object containing the coordinates
    /// @param id the id of the object [default: max of std::size_t]
    explicit Point3dWithID(MathLib::Point3d const& pnt,
        std::size_t id = std::numeric_limits<std::size_t>::max())
        : MathLib::Point3d(pnt), _id(id)
    {}

    /// Default constructor that initializes the id with max of std::size_t
    /// the default constructor of class Point3d.
    Point3dWithID() :
        Point3d(), _id(std::numeric_limits<std::size_t>::max())
    {}

    std::size_t getID() const { return _id; }

protected:
    std::size_t _id;
};

}
