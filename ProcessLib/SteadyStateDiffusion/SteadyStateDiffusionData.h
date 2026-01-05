// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;
};

}  // namespace SteadyStateDiffusion
}  // namespace ProcessLib
