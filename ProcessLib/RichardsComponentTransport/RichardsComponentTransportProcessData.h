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

#include "PorousMediaProperties.h"

namespace MaterialPropertyLib
{
class MaterialSpatialDistributionMap;
}

namespace ProcessLib
{
template <typename ReturnType>
struct Parameter;

namespace RichardsComponentTransport
{
struct RichardsComponentTransportProcessData
{
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
    PorousMediaProperties porous_media_properties;
    Eigen::VectorXd const specific_body_force;
    bool const has_gravity;
};

}  // namespace RichardsComponentTransport
}  // namespace ProcessLib
