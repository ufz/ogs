// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

#include "Base.h"
#include "BaseLib/StrongType.h"

namespace ProcessLib::ThermoRichardsMechanics
{
using LiquidViscosityData =
    BaseLib::StrongType<double, struct LiquidViscosityDataTag>;

constexpr std::string_view ioName(struct LiquidViscosityDataTag*)
{
    return "viscosity";
}
}  // namespace ProcessLib::ThermoRichardsMechanics