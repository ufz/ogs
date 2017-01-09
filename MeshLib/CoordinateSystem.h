/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

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
    explicit CoordinateSystem(const GeoLib::AABB &bbox) : _type(getCoordinateSystem(bbox)) {}

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
    unsigned char getCoordinateSystem(const GeoLib::AABB &bbox) const;

    unsigned char _type;
};

} // MeshLib
