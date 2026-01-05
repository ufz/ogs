// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "MaterialLib/SolidModels/MechanicsBase.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim>
using SolidConstitutiveRelation =
    MaterialLib::Solids::MechanicsBase<DisplacementDim>;
}
}  // namespace ProcessLib::ThermoRichardsMechanics
