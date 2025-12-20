// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include "ConstitutiveSetting.h"
#include "TraitsBase.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim>
struct ConstitutiveTraits
{
    using SolidConstitutiveRelation = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStress_StrainTemperature::SolidConstitutiveRelation<
            DisplacementDim>;

    using StatefulData = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStress_StrainTemperature::StatefulData<DisplacementDim>;

    using StatefulDataPrev = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStress_StrainTemperature::StatefulDataPrev<DisplacementDim>;

    using ConstitutiveData = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStress_StrainTemperature::ConstitutiveData<DisplacementDim>;

    using ConstitutiveTempData = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStress_StrainTemperature::ConstitutiveTempData<
            DisplacementDim>;

    using OutputData = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStress_StrainTemperature::OutputData<DisplacementDim>;

    using ConstitutiveModels = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStress_StrainTemperature::ConstitutiveModels<
            DisplacementDim>;

    template <typename TRMProcessData>
    static ConstitutiveModels createConstitutiveModels(
        TRMProcessData const& process_data,
        SolidConstitutiveRelation const& solid_material)
    {
        return ProcessLib::ThermoRichardsMechanics::
            ConstitutiveStress_StrainTemperature::createConstitutiveModels<
                DisplacementDim>(process_data, solid_material);
    }

    using ConstitutiveSetting = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveStress_StrainTemperature::ConstitutiveSetting<
            DisplacementDim>;
};
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
