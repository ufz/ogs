// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ProcessLib/Reflection/ReflectionData.h"

namespace ProcessLib::ThermoRichardsMechanics
{
struct PorosityData
{
    double phi;

    static auto reflect()
    {
        return ProcessLib::Reflection::reflectWithName("porosity",
                                                       &PorosityData::phi);
    }
};

}  // namespace ProcessLib::ThermoRichardsMechanics
