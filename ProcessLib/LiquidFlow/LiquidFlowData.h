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
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;

    const int gravitational_axis_id;
    const double gravitational_acceleration;
};

}  // namespace LiquidFlow
}  // namespace ProcessLib
