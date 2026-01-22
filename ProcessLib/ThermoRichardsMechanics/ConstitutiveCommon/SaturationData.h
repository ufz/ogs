// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct SaturationDataDeriv
{
    double dS_L_dp_cap;
};

struct SaturationData
{
    double S_L;

    static auto reflect()
    {
        return ProcessLib::Reflection::reflectWithName("saturation",
                                                       &SaturationData::S_L);
    }
};

}  // namespace ProcessLib::ThermoRichardsMechanics
