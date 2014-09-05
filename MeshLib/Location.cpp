/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Location.h"

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

bool operator<(const Location& left, const Location& right)
{
    if (left.mesh_id != right.mesh_id) return left.mesh_id < right.mesh_id;
    if (left.item_type != right.item_type) return left.item_type < right.item_type;
    return left.item_id < right.item_id;
}

std::ostream& operator<<(std::ostream& os, Location const& l)
{
    return os << "(" << l.mesh_id
        << ", " << l.item_type
        << ", " << l.item_id
        << ")";
}

}   // namespace MeshLib
