// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "FreeEnergyDensity.h"
#include "Gravity.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "ProcessLib/ConstitutiveRelations/StressData.h"
#include "SolidDensity.h"
#include "SolidMechanics.h"

namespace ProcessLib::SmallDeformation
{
namespace ConstitutiveRelations
{
/// Data whose state must be tracked by the process.
template <int DisplacementDim>
using StatefulData = std::tuple<StressData<DisplacementDim>>;

template <int DisplacementDim>
using StatefulDataPrev = ProcessLib::ConstitutiveRelations::PrevStateOf<
    StatefulData<DisplacementDim>>;

/// Data that is needed for output purposes, but not directly for the assembly.
template <int DisplacementDim>
using OutputData =
    std::tuple<StrainData<DisplacementDim>, FreeEnergyDensityData>;

/// Data that is needed for the equation system assembly.
template <int DisplacementDim>
using ConstitutiveData =
    std::tuple<SolidMechanicsDataStateless<DisplacementDim>,
               VolumetricBodyForce<DisplacementDim>>;

/// Data that stores intermediate values, which are not needed outside the
/// constitutive setting.
template <int DisplacementDim>
using ConstitutiveTempData =
    std::tuple<PrevState<StrainData<DisplacementDim>>, SolidDensity>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::SmallDeformation
