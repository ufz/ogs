/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TEquation.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
void FT1Model::eval(
    double const dt,
    InternalEnergyData const& internal_energy_data,
    PrevState<InternalEnergyData> const& internal_energy_data_prev,
    FT1Data& fT_1) const
{
    if (dt == 0)
    {
        fT_1.m = 0;
        return;
    }

    auto const rho_u_eff_dot =
        (internal_energy_data() - **internal_energy_data_prev) / dt;
    fT_1.m = rho_u_eff_dot;
}

void FT1Model::dEval(double const dt,
                     EffectiveVolumetricInternalEnergyDerivatives const&
                         effective_volumetric_internal_energy_d_data,
                     FT1DerivativeData& dfT_1) const
{
    if (dt == 0)
    {
        dfT_1.dp_GR = 0;
        dfT_1.dp_cap = 0;
        dfT_1.dT = 0;
        return;
    }

    dfT_1.dp_GR =
        effective_volumetric_internal_energy_d_data.drho_u_eff_dp_GR / dt;

    dfT_1.dp_cap =
        effective_volumetric_internal_energy_d_data.drho_u_eff_dp_cap / dt;

    dfT_1.dT = effective_volumetric_internal_energy_d_data.drho_u_eff_dT / dt;
}

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
