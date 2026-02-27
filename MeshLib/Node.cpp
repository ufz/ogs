// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MeshLib/Node.h"

#include "Elements/Element.h"

namespace MeshLib
{
Node::Node(const double coords[3], std::size_t id)
    : MathLib::Point3dWithID(
          std::array<double, 3>{{coords[0], coords[1], coords[2]}}, id)
{
}

Node::Node(std::array<double, 3> const& coords, std::size_t id)
    : MathLib::Point3dWithID(coords, id)
{
}

Node::Node(double x, double y, double z, std::size_t id)
    : MathLib::Point3dWithID(std::array<double, 3>({{x, y, z}}), id)
{
}

std::ostream& operator<<(std::ostream& os, MeshLib::Node const& n)
{
    auto const default_precision = os.precision();
    os << "node #" << n.getID() << " { "
       << std::setprecision(std::numeric_limits<double>::digits10 + 1)
       << static_cast<MathLib::Point3d const&>(n) << "}";
    os << std::setprecision(default_precision);
    return os;
}
}  // namespace MeshLib
