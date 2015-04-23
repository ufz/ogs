/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef COORDINATESYSTEMTYPE_H_
#define COORDINATESYSTEMTYPE_H_

#include <cassert>
#include <cstddef>
#include <cmath>

#include "GeoLib/AABB.h"

namespace MeshLib
{

class Element;

/**
 * \brief Coordinate system type
 */
enum class CoordinateSystemType
{
    X = 10,
    Y = 11,
    Z = 12,
    XY = 21,
    XZ = 22,
    YZ = 23,
    XYZ = 32
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
    explicit CoordinateSystem(CoordinateSystemType coord) : _type (coord) {}

    /// Decides for a given element
    explicit CoordinateSystem(const Element &ele);

    /// Decides a coordinate system from a bounding box
    template <class T>
    explicit CoordinateSystem(const GeoLib::AABB<T> &bbox) : _type(getCoordinateSystem(bbox)) {}

    /// get this coordinate type
    CoordinateSystemType getType() const {
        return _type;
    }

    /// get dimension size
    unsigned getDimension() const {
        switch (_type) {
            case CoordinateSystemType::X:
                return 1;
            case CoordinateSystemType::Y:
            case CoordinateSystemType::XY:
                return 2;
            default:
                return 3;
        }
    }

    /// has X dimension
    bool hasX() const;

    /// has Y dimension
    bool hasY() const;

    /// has z dimension
    bool hasZ() const;

private:
    template <class T>
    CoordinateSystemType getCoordinateSystem(const GeoLib::AABB<T> &bbox) const;

    CoordinateSystemType _type;
};

template <class T>
CoordinateSystemType CoordinateSystem::getCoordinateSystem(const GeoLib::AABB<T> &bbox) const
{
    MeshLib::CoordinateSystemType coords;

    T pt_diff = bbox.getMaxPoint() - bbox.getMinPoint();
    bool hasX = std::abs(pt_diff[0]) > .0;
    bool hasY = std::abs(pt_diff[1]) > .0;
    bool hasZ = std::abs(pt_diff[2]) > .0;

    if (hasX) {
        if (hasY) {
            if (hasZ) {
                coords = CoordinateSystemType::XYZ;
            } else {
                coords = CoordinateSystemType::XY;
            }
        } else if (hasZ) {
            coords = CoordinateSystemType::XZ;
        } else {
            coords = CoordinateSystemType::X;
        }
    } else if (hasY) {
        if (hasZ) {
            coords = CoordinateSystemType::YZ;
        } else {
            coords = CoordinateSystemType::Y;
        }
    } else if (hasZ) {
        coords = CoordinateSystemType::Z;
    } else {
        coords = CoordinateSystemType::X;
    }
    return coords;
}

} // MeshLib

#endif // COORDINATESYSTEMTYPE_H_
