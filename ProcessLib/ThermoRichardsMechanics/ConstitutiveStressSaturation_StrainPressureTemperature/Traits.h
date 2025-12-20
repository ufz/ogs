// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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

    template <typename TRMProcessData>
    static ConstitutiveModels createConstitutiveModels(
        TRMProcessData const& process_data,
        SolidConstitutiveRelation const& solid_material)
    {
        return ProcessLib::ThermoRichardsMechanics::
            ConstitutiveStressSaturation_StrainPressureTemperature::
                createConstitutiveModels<DisplacementDim>(process_data,
                                                          solid_material);
    }

    using ConstitutiveSetting = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStressSaturation_StrainPressureTemperature::
            ConstitutiveSetting<DisplacementDim>;
};
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
