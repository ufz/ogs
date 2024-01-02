/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
