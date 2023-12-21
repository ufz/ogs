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
#include "ProcessLib/Reflection/ReflectionData.h"
#include "SolidDensity.h"
#include "SolidMechanics.h"

namespace ProcessLib::LargeDeformation
{
namespace ConstitutiveRelations
{
/// Data whose state must be tracked by the process.
template <int DisplacementDim>
struct StatefulData
{
    StressData<DisplacementDim> stress_data;

    static auto reflect()
    {
        using Self = StatefulData<DisplacementDim>;

        return Reflection::reflectWithoutName(&Self::stress_data);
    }
};

/// Data whose state must be tracked by the process.
template <int DisplacementDim>
struct StatefulDataPrev
{
    PrevState<StressData<DisplacementDim>> stress_data;

    StatefulDataPrev<DisplacementDim>& operator=(
        StatefulData<DisplacementDim> const& state)
    {
        stress_data = state.stress_data;

        return *this;
    }
};

/// Data that is needed for output purposes, but not directly for the assembly.
template <int DisplacementDim>
struct OutputData
{
    StrainData<DisplacementDim> eps_data;
    DeformationGradientData<DisplacementDim> deformation_gradient_data;

    static auto reflect()
    {
        using Self = OutputData<DisplacementDim>;

        return Reflection::reflectWithoutName(&Self::eps_data,
                                              &Self::deformation_gradient_data);
    }
};

/// Data that is needed for the equation system assembly.
template <int DisplacementDim>
struct ConstitutiveData
{
    SolidMechanicsDataStateless<DisplacementDim> s_mech_data;
    VolumetricBodyForce<DisplacementDim> volumetric_body_force;
};

/// Data that stores intermediate values, which are not needed outside the
/// constitutive setting.
template <int DisplacementDim>
struct ConstitutiveTempData
{
    PrevState<StrainData<DisplacementDim>> eps_data_prev;
    PrevState<DeformationGradientData<DisplacementDim>>
        deformation_gradient_data_prev;
    SolidDensity rho_SR;
};
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::LargeDeformation
