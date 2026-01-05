// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

#include "BaseLib/StrongType.h"

namespace ProcessLib::RichardsMechanics
{
// Apparent dry solid density
using DrySolidDensity = BaseLib::StrongType<double, struct DrySolidDensityTag>;

constexpr std::string_view ioName(struct DrySolidDensityTag*)
{
    return "dry_density_solid";
}
}  // namespace ProcessLib::RichardsMechanics
