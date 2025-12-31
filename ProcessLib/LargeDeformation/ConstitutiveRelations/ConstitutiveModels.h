// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Gravity.h"
#include "ProcessLib/Graph/ConstructModels.h"
#include "SolidDensity.h"
#include "SolidMechanics.h"
#include "SpecificBodyForceData.h"

namespace ProcessLib::LargeDeformation
{
namespace ConstitutiveRelations
{
/// Constitutive models used for assembly.
template <int DisplacementDim>
using ConstitutiveModels = std::tuple<SolidMechanicsModel<DisplacementDim>,
                                      SolidDensityModel,
                                      GravityModel<DisplacementDim>>;

template <int DisplacementDim, typename LDProcessData>
ConstitutiveModels<DisplacementDim> createConstitutiveModels(
    LDProcessData const& process_data,
    SolidConstitutiveRelation<DisplacementDim> const& solid_material)
{
    return ProcessLib::Graph::constructModels<
        ConstitutiveModels<DisplacementDim>>(
        SpecificBodyForceData<DisplacementDim>{
            process_data.specific_body_force},
        solid_material);
}
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::LargeDeformation
