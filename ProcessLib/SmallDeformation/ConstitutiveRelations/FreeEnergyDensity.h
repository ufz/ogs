// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::SmallDeformation
{
struct FreeEnergyDensityData
{
    double free_energy_density;

    static auto reflect()
    {
        return ProcessLib::Reflection::reflectWithName(
            "free_energy_density", &FreeEnergyDensityData::free_energy_density);
    }
};
}  // namespace ProcessLib::SmallDeformation
