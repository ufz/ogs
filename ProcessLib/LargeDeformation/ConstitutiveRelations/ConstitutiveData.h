// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Gravity.h"
#include "ProcessLib/ConstitutiveRelations/Base.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "ProcessLib/ConstitutiveRelations/StressData.h"
#include "SolidDensity.h"
#include "SolidMechanics.h"

namespace ProcessLib::LargeDeformation
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
/// This tuple is populated by the FEM code after constitutive evaluation.
template <int DisplacementDim>
using OutputData = std::tuple<StrainData<DisplacementDim>,
                              DeformationGradientData<DisplacementDim>>;

/// Data that is needed for the equation system assembly.
template <int DisplacementDim>
using ConstitutiveData =
    std::tuple<SolidMechanicsDataStateless<DisplacementDim>,
               VolumetricBodyForce<DisplacementDim>>;

/// Data that stores intermediate values, which are not needed outside the
/// constitutive setting.
template <int DisplacementDim>
using ConstitutiveTempData =
    std::tuple<PrevState<DeformationGradientData<DisplacementDim>>,
               SolidDensity>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::LargeDeformation
