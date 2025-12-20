// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"

namespace ProcessLib::HeatConduction
{
struct HeatConductionProcessData
{
    MaterialPropertyLib::MaterialSpatialDistributionMap media_map;

    /// If set mass lumping will be applied to the equation.
    bool const mass_lumping;

    int const mesh_space_dimension;
};
}  // namespace ProcessLib::HeatConduction
