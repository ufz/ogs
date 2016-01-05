/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef COORDINATESYSTEMTYPE_H_
#define COORDINATESYSTEMTYPE_H_

#include <cmath>

#include "GeoLib/AABB.h"
#include "MathLib/Vector3.h"

namespace MeshLib
{

class Element;

/**
 * \brief Coordinate system type
 */
struct CoordinateSystemType
{
    enum type
    {
        X = 0x01,
        Y = 0x02,
        Z = 0x04
    };
};

/**
 * \brief Coordinate systems
 *
 *
 */
class CoordinateSystem
{
public:
    /// User provided coordinate system
    explicit CoordinateSystem(unsigned char coord) : _type (coord) {}

    /// Decides for a given element
    explicit CoordinateSystem(const Element &ele);

    /// Decides a coordinate system from a bounding box
    template <class T>
    explicit CoordinateSystem(const GeoLib::AABB<T> &bbox) : _type(getCoordinateSystem(bbox)) {}

    /// get this coordinate type
    unsigned char getType() const { return _type; }

    /// get dimension size
    unsigned getDimension() const {
        if (hasZ())
            return 3;
        else if (hasY())
            return 2;
        else
            return 1;
    }

    /// has X dimension
    bool hasX() const { return (_type & CoordinateSystemType::type::X) != 0; }

    /// has Y dimension
    bool hasY() const { return (_type & CoordinateSystemType::type::Y) != 0; }

    /// has z dimension
    bool hasZ() const { return (_type & CoordinateSystemType::type::Z) != 0; }

private:
    template <class T>
    unsigned char getCoordinateSystem(const GeoLib::AABB<T> &bbox) const;

    unsigned char _type;
};

template <class T>
unsigned char CoordinateSystem::getCoordinateSystem(const GeoLib::AABB<T> &bbox) const
{
    unsigned char coords = 0;

    const MathLib::Vector3 pt_diff(bbox.getMinPoint(), bbox.getMaxPoint());

    // The axis aligned bounding box is a from the right half-open interval.
    // Therefore, the difference between the particular coordinates of the
    // points is modified by the unit in the last place towards zero.
    if (std::nexttoward(std::abs(pt_diff[0]), 0.0) > .0)
        coords |= CoordinateSystemType::X;
    if (std::nexttoward(std::abs(pt_diff[1]), 0.0) > .0)
        coords |= CoordinateSystemType::Y;
    if (std::nexttoward(std::abs(pt_diff[2]), 0.0) > .0)
        coords |= CoordinateSystemType::Z;

    return coords;
}

} // MeshLib

#endif // COORDINATESYSTEMTYPE_H_
