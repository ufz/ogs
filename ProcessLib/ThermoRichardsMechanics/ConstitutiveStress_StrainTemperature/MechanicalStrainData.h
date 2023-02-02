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

#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Base.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStress_StrainTemperature
{
template <int DisplacementDim>
struct MechanicalStrainData
{
    // TODO it seems fragile that some data have to be initialized that way.
    KelvinVector<DisplacementDim> eps_m = KV::KVzero<DisplacementDim>();

    static auto reflect()
    {
        using Self = MechanicalStrainData<DisplacementDim>;

        return ProcessLib::Reflection::reflectWithName("eps_m", &Self::eps_m);
    }
};
}  // namespace ConstitutiveStress_StrainTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
