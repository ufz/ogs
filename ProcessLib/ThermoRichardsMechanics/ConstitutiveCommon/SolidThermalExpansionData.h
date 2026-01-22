// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
struct SolidThermalExpansionData
{
    KelvinVector<DisplacementDim> solid_linear_thermal_expansivity_vector;
};

}  // namespace ProcessLib::ThermoRichardsMechanics
