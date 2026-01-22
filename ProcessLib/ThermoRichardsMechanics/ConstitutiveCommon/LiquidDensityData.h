// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct LiquidDensityData
{
    double rho_LR;
    double drho_LR_dp;
    double drho_LR_dT;

    static auto reflect()
    {
        return ::ProcessLib::Reflection::reflectWithName(
            "liquid_density", &LiquidDensityData::rho_LR);
    }
};

}  // namespace ProcessLib::ThermoRichardsMechanics
