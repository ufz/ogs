// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MechanicalStrainModel.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
void MechanicalStrainModel<DisplacementDim>::eval(
    TemperatureData const& T_data,
    SolidThermalExpansionData<DisplacementDim> const& s_therm_exp_data,
    ProcessLib::ConstitutiveRelations::StrainData<DisplacementDim> const&
        strain_data,
    KelvinVector<DisplacementDim> const& eps_prev,
    PrevState<MechanicalStrainData<DisplacementDim>> const& eps_m_prev,
    SwellingDataStateless<DisplacementDim> const& swelling_data,
    MechanicalStrainData<DisplacementDim>& out) const
{
    MathLib::KelvinVector::KelvinVectorType<DisplacementDim> const
        dthermal_strain =
            s_therm_exp_data.solid_linear_thermal_expansivity_vector *
            (T_data.T - T_data.T_prev);

    out.eps_m.noalias() = eps_m_prev->eps_m + strain_data.eps - eps_prev -
                          dthermal_strain + swelling_data.eps_m;
}

template struct MechanicalStrainModel<2>;
template struct MechanicalStrainModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
