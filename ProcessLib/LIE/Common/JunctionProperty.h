// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>
#include <array>

#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace LIE
{
struct JunctionProperty final
{
    JunctionProperty(int const junction_id_,
                     MeshLib::Node const& junctionNode,
                     std::array<int, 2> const fracture_ids_)
        : coords{junctionNode.data()},
          node_id{junctionNode.getID()},
          fracture_ids{fracture_ids_},
          junction_id{junction_id_}
    {
    }
    Eigen::Vector3d const coords;
    std::size_t const node_id;
    std::array<int, 2> const fracture_ids;
    int const junction_id;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace LIE
}  // namespace ProcessLib
