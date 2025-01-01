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
// Apparent dry solid density
using DrySolidDensity = BaseLib::StrongType<double, struct DrySolidDensityTag>;

constexpr std::string_view ioName(struct DrySolidDensityTag*)
{
    return "dry_density_solid";
}
}  // namespace ProcessLib::RichardsMechanics
