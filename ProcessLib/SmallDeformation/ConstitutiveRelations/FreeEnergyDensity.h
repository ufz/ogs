/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

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
