/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <array>

#include <Eigen/Eigen>

namespace ProcessLib
{
namespace LIE
{
struct JunctionProperty final
{
    Eigen::Vector3d coords;
    int junction_id;
    int node_id;
    std::array<int, 2> fracture_IDs;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace LIE
}  // namespace ProcessLib
