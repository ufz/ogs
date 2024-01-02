/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Gravity.h"
#include "SolidDensity.h"
#include "SolidMechanics.h"

namespace ProcessLib::SmallDeformation
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
}  // namespace ProcessLib::SmallDeformation
