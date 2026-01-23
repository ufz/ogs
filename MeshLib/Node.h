// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <spdlog/fmt/ostr.h>

#include "MathLib/Point3dWithID.h"

namespace MeshToolsLib
{
class MeshRevision;
}

namespace MeshLib
{
/**
 * A mesh node with coordinates in 3D space.
 */
class Node : public MathLib::Point3dWithID
{
    /* friend classes: */
    friend class Mesh;
    friend class MeshToolsLib::MeshRevision;

public:
    /// Constructor using a coordinate array
    explicit Node(const double coords[3],
                  std::size_t id = std::numeric_limits<std::size_t>::max());

    /// Constructor using a coordinate array
    explicit Node(std::array<double, 3> const& coords,
                  std::size_t id = std::numeric_limits<std::size_t>::max());

    /// Constructor using single coordinates
    Node(double x,
         double y,
         double z,
         std::size_t id = std::numeric_limits<std::size_t>::max());

    friend std::ostream& operator<<(std::ostream& os, Node const& n);
}; /* class */
}  // namespace MeshLib

namespace fmt
{
template <>
struct formatter<::MeshLib::Node> : ostream_formatter
{
};
}  // namespace fmt
