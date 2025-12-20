// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/StrongType.h"

namespace ProcessLib::RichardsMechanics
{
// effective density of both the solid and fluid phases
using Density = BaseLib::StrongType<double, struct DensityTag>;
}  // namespace ProcessLib::RichardsMechanics
