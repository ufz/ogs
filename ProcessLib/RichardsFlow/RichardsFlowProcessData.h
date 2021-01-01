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

#include <memory>

namespace MaterialPropertyLib
{
class MaterialSpatialDistributionMap;
}

namespace ProcessLib
{
namespace RichardsFlow
{
struct RichardsFlowProcessData
{
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
    bool const has_mass_lumping;
};

}  // namespace RichardsFlow
}  // namespace ProcessLib
