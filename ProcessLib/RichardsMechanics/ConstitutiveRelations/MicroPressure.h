// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

#include "BaseLib/StrongType.h"

namespace ProcessLib::RichardsMechanics
{
using MicroPressure = BaseLib::StrongType<double, struct MicroPressureTag>;

constexpr std::string_view ioName(struct MicroPressureTag*)
{
    return "micro_pressure";
}
}  // namespace ProcessLib::RichardsMechanics
