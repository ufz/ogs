// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <string_view>

#include "BaseLib/StrongType.h"

namespace ProcessLib::RichardsMechanics
{
using MicroSaturation = BaseLib::StrongType<double, struct MicroSaturationTag>;

constexpr std::string_view ioName(struct MicroSaturationTag*)
{
    return "micro_saturation";
}
}  // namespace ProcessLib::RichardsMechanics
