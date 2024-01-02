/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Core>

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"

namespace ProcessLib
{
namespace RichardsComponentTransport
{
struct RichardsComponentTransportProcessData
{
    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
};

}  // namespace RichardsComponentTransport
}  // namespace ProcessLib
