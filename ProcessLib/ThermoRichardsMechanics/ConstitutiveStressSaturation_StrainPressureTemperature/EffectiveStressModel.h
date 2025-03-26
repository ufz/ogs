/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "ProcessLib/ConstitutiveRelations/StressData.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Biot.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Bishops.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TotalStressData.h"
#include "TraitsBase.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
template <int DisplacementDim>
struct EffectiveStressModel
{
    void eval(
        CapillaryPressureData<DisplacementDim> const& p_cap_data,
        BiotData const& biot_data,
        BishopsData const& bishops_data,
        TotalStressData<DisplacementDim> const& total_stress_data,
        ProcessLib::ConstitutiveRelations::EffectiveStressData<DisplacementDim>&
            sigma_eff_data) const;
};

extern template struct EffectiveStressModel<2>;
extern template struct EffectiveStressModel<3>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
