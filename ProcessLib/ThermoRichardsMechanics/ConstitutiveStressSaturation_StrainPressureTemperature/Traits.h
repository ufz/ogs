/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ConstitutiveSetting.h"
#include "TraitsBase.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
template <int DisplacementDim>
struct ConstitutiveTraits
{
    using SolidConstitutiveRelation = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStressSaturation_StrainPressureTemperature::
            SolidConstitutiveRelation<DisplacementDim>;

    using StatefulData = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStressSaturation_StrainPressureTemperature::StatefulData<
            DisplacementDim>;

    using StatefulDataPrev = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStressSaturation_StrainPressureTemperature::
            StatefulDataPrev<DisplacementDim>;

    using ConstitutiveData = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStressSaturation_StrainPressureTemperature::
            ConstitutiveData<DisplacementDim>;

    using ConstitutiveTempData = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStressSaturation_StrainPressureTemperature::
            ConstitutiveTempData<DisplacementDim>;

    using OutputData = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStressSaturation_StrainPressureTemperature::OutputData<
            DisplacementDim>;

    using ConstitutiveModels = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStressSaturation_StrainPressureTemperature::
            ConstitutiveModels<DisplacementDim>;

    using ConstitutiveSetting = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStressSaturation_StrainPressureTemperature::
            ConstitutiveSetting<DisplacementDim>;
};
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
