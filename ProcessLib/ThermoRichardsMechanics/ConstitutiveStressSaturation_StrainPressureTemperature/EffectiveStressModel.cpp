// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "EffectiveStressModel.h"

namespace ProcessLib::ThermoRichardsMechanics
{
namespace ConstitutiveStressSaturation_StrainPressureTemperature
{
template <int DisplacementDim>
void EffectiveStressModel<DisplacementDim>::eval(
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    BiotData const& biot_data,
    BishopsData const& bishops_data,
    TotalStressData<DisplacementDim> const& total_stress_data,
    ProcessLib::ConstitutiveRelations::EffectiveStressData<DisplacementDim>&
        sigma_eff_data) const
{
    auto const& identity2 = MathLib::KelvinVector::Invariants<
        MathLib::KelvinVector::kelvin_vector_dimensions(
            DisplacementDim)>::identity2;

    sigma_eff_data.sigma_eff.noalias() =
        total_stress_data.sigma_total -
        biot_data() * bishops_data.chi_S_L * p_cap_data.p_cap * identity2;
}

template struct EffectiveStressModel<2>;
template struct EffectiveStressModel<3>;
}  // namespace ConstitutiveStressSaturation_StrainPressureTemperature
}  // namespace ProcessLib::ThermoRichardsMechanics
