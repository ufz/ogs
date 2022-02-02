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

#include <memory>
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"

namespace ProcessLib::HeatConduction
{
struct HeatConductionProcessData
{
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;

    /// If set mass lumping will be applied to the equation.
    bool const mass_lumping;

    int const mesh_space_dimension;
};
}  // namespace ProcessLib::HeatConduction
