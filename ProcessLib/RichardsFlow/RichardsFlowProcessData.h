// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Eigen/Core>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"

namespace ProcessLib
{
namespace RichardsFlow
{
struct RichardsFlowProcessData
{
    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
    bool const has_mass_lumping;
};

}  // namespace RichardsFlow
}  // namespace ProcessLib
