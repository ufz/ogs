/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ConstitutiveSetting.h"
#include "TraitsBase.h"

namespace ProcessLib::ThermoRichardsMechanics::ConstitutiveOriginal
{
template <int DisplacementDim>
struct ConstitutiveTraits
{
    using SolidConstitutiveRelation = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveOriginal::SolidConstitutiveRelation<DisplacementDim>;

    using StatefulData =
        ProcessLib::ThermoRichardsMechanics::ConstitutiveOriginal::StatefulData<
            DisplacementDim>;

    using ConstitutiveData = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveOriginal::ConstitutiveData<DisplacementDim>;

    using ConstitutiveTempData = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveOriginal::ConstitutiveTempData<DisplacementDim>;

    using OutputData =
        ProcessLib::ThermoRichardsMechanics::ConstitutiveOriginal::OutputData<
            DisplacementDim>;

    using ConstitutiveModels = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveOriginal::ConstitutiveModels<DisplacementDim>;

    using ConstitutiveSetting = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveOriginal::ConstitutiveSetting<DisplacementDim>;

    using ElasticTangentStiffnessModel = ProcessLib::ThermoRichardsMechanics::
        ConstitutiveOriginal::ElasticTangentStiffnessModel<DisplacementDim>;
};
}  // namespace ProcessLib::ThermoRichardsMechanics::ConstitutiveOriginal
