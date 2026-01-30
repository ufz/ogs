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
// Explicit instantiation declarations to avoid multiple-definition issues.
extern template struct SolidThermalExpansionData<2>;
extern template struct SolidThermalExpansionData<3>;

}  // namespace ProcessLib::ThermoRichardsMechanics
