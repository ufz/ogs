/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Location.h"
#include <ostream>

namespace MeshLib
{

std::ostream& operator<<(std::ostream& os, MeshItemType const& t)
{
    switch (t)
    {
    case MeshItemType::Node: return os << "N";
    case MeshItemType::Edge: return os << "E";
    case MeshItemType::Face: return os << "F";
    case MeshItemType::Cell: return os << "C";
    };
    return os;
}

std::ostream& operator<<(std::ostream& os, Location const& l)
{
    return os << "(" << l.mesh_id
        << ", " << l.item_type
        << ", " << l.item_id
        << ")";
}

}   // namespace MeshLib
