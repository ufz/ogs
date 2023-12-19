/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Biot.h"
#include "ElasticTangentStiffnessModel.h"
#include "SolidCompressibility.h"
#include "SolidMechanics.h"
#include "SolidThermalExpansion.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
/// Constitutive models used for assembly.
template <int DisplacementDim>
struct ConstitutiveModels
{
    explicit ConstitutiveModels(
        SolidConstitutiveRelation<DisplacementDim> const& solid_material)
        : elastic_tangent_stiffness_model(solid_material),
          beta_p_SR_model(solid_material)
    {
    }

    ElasticTangentStiffnessModel<DisplacementDim>
        elastic_tangent_stiffness_model;
    BiotModel biot_model;
    SolidCompressibilityModel<DisplacementDim,
                              SolidConstitutiveRelation<DisplacementDim>>
        beta_p_SR_model;
    SolidThermalExpansionModel<DisplacementDim> s_therm_exp_model;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
