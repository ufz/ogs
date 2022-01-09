/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>
#include <memory>

namespace MaterialPropertyLib
{
class MaterialSpatialDistributionMap;
}

namespace ProcessLib
{
namespace LiquidFlow
{
struct LiquidFlowData final
{
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;

    /// A vector of the rotation matrices for all elements.
    std::vector<Eigen::MatrixXd> const element_rotation_matrices;

    int const mesh_space_dimension;

    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
};

}  // namespace LiquidFlow
}  // namespace ProcessLib
