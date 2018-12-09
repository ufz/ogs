/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Eigen>

#include "MeshLib/Node.h"

namespace ProcessLib
{
namespace LIE
{
struct BranchProperty final
{
    BranchProperty(MeshLib::Node const& branchNode,
                   int const master_fracture_id_,
                   int const slave_fracture_id_)
        : coords{branchNode.getCoords()},
          node_id{branchNode.getID()},
          master_fracture_id{master_fracture_id_},
          slave_fracture_id{slave_fracture_id_}
    {
    }

    Eigen::Vector3d const coords;
    // unit vector normal to the master fracture in a direction to the slave
    Eigen::Vector3d normal_vector_branch;
    std::size_t const node_id;
    int const master_fracture_id;
    int const slave_fracture_id;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace LIE
}  // namespace ProcessLib
