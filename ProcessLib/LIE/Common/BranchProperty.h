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

namespace ProcessLib
{
namespace LIE
{
struct BranchProperty final
{
    Eigen::Vector3d coords;
    // unit vector normal to the master fracture in a direction to the slave
    Eigen::Vector3d normal_vector_branch;
    int node_id;
    int master_fracture_ID;
    int slave_fracture_ID;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace LIE
}  // namespace ProcessLib
