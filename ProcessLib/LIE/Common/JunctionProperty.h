/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>
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
        : coords{junctionNode.getCoords()},
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
