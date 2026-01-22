// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct SolidDensityData
{
    double rho_SR;
    double dry_density_solid;

    static auto reflect()
    {
        return ProcessLib::Reflection::reflectWithName(
            "dry_density_solid", &SolidDensityData::dry_density_solid);
    }
};

}  // namespace ProcessLib::ThermoRichardsMechanics
