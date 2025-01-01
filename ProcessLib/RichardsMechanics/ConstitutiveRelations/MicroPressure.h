/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
