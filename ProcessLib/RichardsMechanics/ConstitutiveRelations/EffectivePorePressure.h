// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "BaseLib/StrongType.h"

namespace ProcessLib::RichardsMechanics
{
using EffectivePorePressure =
    BaseLib::StrongType<double, struct EffectivePorePressureTag>;
}  // namespace ProcessLib::RichardsMechanics
