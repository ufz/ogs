// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "Gravity.h"
#include "SolidDensity.h"
#include "SolidMechanics.h"

namespace ProcessLib::LargeDeformation
{
namespace ConstitutiveRelations
{
/// Constitutive models used for assembly.
template <int DisplacementDim>
struct ConstitutiveModels
{
    template <typename TRMProcessData>
    explicit ConstitutiveModels(
        TRMProcessData const& process_data,
        SolidConstitutiveRelation<DisplacementDim> const& solid_material)
        : s_mech_model(solid_material),
          gravity_model(process_data.specific_body_force)
    {
    }

    SolidMechanicsModel<DisplacementDim> s_mech_model;
    SolidDensityModel rho_S_model;
    GravityModel<DisplacementDim> gravity_model;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::LargeDeformation
