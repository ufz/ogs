// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "ThermoOsmosisModel.h"

#include "ProcessLib/Common/ThermoOsmosis/ThermoOsmoticCoefficient.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void ThermoOsmosisModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    TemperatureData<DisplacementDim> const& T_data,
    LiquidDensityData const& rho_L_data,
    PermeabilityData<DisplacementDim> const& perm_data,
    LiquidViscosityData const& mu_L_data,
    ThermoOsmosisData<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;
    // Holds no primary variables: the properties read below are evaluated by
    // the models this one takes its data from, not here.
    MPL::VariableArray variables;

    // Ki is row-major, the helper takes the default column-major layout.
    Eigen::Matrix<double, DisplacementDim, DisplacementDim> const
        intrinsic_permeability = perm_data.Ki;

    Eigen::Matrix<double, DisplacementDim, DisplacementDim> const
        K_pT_thermal_osmosis =
            ProcessLib::getThermoOsmoticCoefficient<DisplacementDim>(
                media_data.medium, variables, x_t.x, x_t.t, x_t.dt,
                intrinsic_permeability, *mu_L_data);

    out.K_pT_Laplace = rho_L_data.rho_LR * K_pT_thermal_osmosis;

    out.K_Tp_Laplace = T_data.T * K_pT_thermal_osmosis;

    out.seepage_velocity_contribution = -K_pT_thermal_osmosis * T_data.grad_T;
}

template struct ThermoOsmosisModel<2>;
template struct ThermoOsmosisModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
