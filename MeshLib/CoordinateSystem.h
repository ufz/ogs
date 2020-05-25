/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    explicit CoordinateSystem(unsigned char coord) : type_ (coord) {}

    /// Decides for a given element
    explicit CoordinateSystem(const Element &ele);

    /// Decides a coordinate system from a bounding box
    explicit CoordinateSystem(const GeoLib::AABB &bbox) : type_(getCoordinateSystem(bbox)) {}

    /// get this coordinate type
    unsigned char getType() const { return type_; }

    /// get dimension size
    unsigned getDimension() const {
        if (hasZ())
        {
            return 3;
        }
        if (hasY())
        {
            return 2;
        }

        return 1;
    }

    /// has X dimension
    bool hasX() const { return (type_ & CoordinateSystemType::type::X) != 0; }

    /// has Y dimension
    bool hasY() const { return (type_ & CoordinateSystemType::type::Y) != 0; }

    /// has z dimension
    bool hasZ() const { return (type_ & CoordinateSystemType::type::Z) != 0; }

private:
    unsigned char getCoordinateSystem(const GeoLib::AABB &bbox) const;

    unsigned char type_;
};

}  // namespace MeshLib
