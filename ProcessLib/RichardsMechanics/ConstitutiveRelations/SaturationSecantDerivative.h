// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Base.h"

namespace ProcessLib::RichardsMechanics
{
struct SaturationSecantDerivative
{
    double DeltaS_L_Deltap_cap = nan;
};
}  // namespace ProcessLib::RichardsMechanics
