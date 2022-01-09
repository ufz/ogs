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

namespace ProcessLib
{

template <typename T>
struct Parameter;

namespace SteadyStateDiffusion
{
struct SteadyStateDiffusionData final
{
    std::unique_ptr<MaterialPropertyLib::MaterialSpatialDistributionMap>
        media_map;
};

} // namespace SteadyStateDiffusion
} // namespace ProcessLib
