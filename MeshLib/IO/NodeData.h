// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <cstddef>

namespace MeshLib::IO
{
/// struct NodeData used for parallel reading and also partitioning
struct NodeData
{
    NodeData() = default;

    NodeData(std::size_t id, double ix, double iy, double iz)
        : index(id), x(ix), y(iy), z(iz)
    {
    }

    std::size_t index;  ///< Global node index.
    double x;
    double y;
    double z;
};
}  // end namespace MeshLib::IO
