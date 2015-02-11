/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CoordinateSystem.h"

#include <limits>

namespace MeshLib
{

bool CoordinateSystem::hasX() const {
    switch (_type) {
    case CoordinateSystemType::X:
    case CoordinateSystemType::XY:
    case CoordinateSystemType::XZ:
    case CoordinateSystemType::XYZ:
        return true;
    default:
        return false;
    }
}

bool CoordinateSystem::hasY() const {
    switch (_type) {
    case CoordinateSystemType::Y:
    case CoordinateSystemType::XY:
    case CoordinateSystemType::YZ:
    case CoordinateSystemType::XYZ:
        return true;
    default:
        return false;
    }
}

bool CoordinateSystem::hasZ() const {
    switch (_type) {
    case CoordinateSystemType::Z:
    case CoordinateSystemType::XZ:
    case CoordinateSystemType::YZ:
    case CoordinateSystemType::XYZ:
        return true;
    default:
        return false;
    }
}


} // end
