// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "Location.h"

#include <ostream>

namespace MeshLib
{
std::ostream& operator<<(std::ostream& os, MeshItemType const& t)
{
    switch (t)
    {
        case MeshItemType::Node:
            return os << "N";
        case MeshItemType::Edge:
            return os << "E";
        case MeshItemType::Face:
            return os << "F";
        case MeshItemType::Cell:
            return os << "C";
        case MeshItemType::IntegrationPoint:
            return os << "I";
    };
    return os;
}

std::ostream& operator<<(std::ostream& os, Location const& l)
{
    return os << "(" << l.mesh_id << ", " << l.item_type << ", " << l.item_id
              << ")";
}

}  // namespace MeshLib
