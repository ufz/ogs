/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ThermoOsmosis.h"

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void ThermoOsmosisModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    TemperatureData<DisplacementDim> const& T_data,
    LiquidDensityData const& rho_L_data,
    ThermoOsmosisData<DisplacementDim>& out) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;

    auto const& solid_phase = media_data.solid;

    auto const K_pT_thermal_osmosis =
        (solid_phase.hasProperty(
             MaterialPropertyLib::PropertyType::thermal_osmosis_coefficient)
             ? MaterialPropertyLib::formEigenTensor<DisplacementDim>(
                   solid_phase[MPL::PropertyType::thermal_osmosis_coefficient]
                       .value(variables, x_t.x, x_t.t, x_t.dt))
             : Eigen::MatrixXd::Zero(DisplacementDim, DisplacementDim));

    out.K_pT_Laplace = rho_L_data.rho_LR * K_pT_thermal_osmosis;

    out.K_Tp_Laplace = T_data.T * K_pT_thermal_osmosis;

    out.seepage_velocity_contribution = -K_pT_thermal_osmosis * T_data.grad_T;
}

template struct ThermoOsmosisModel<2>;
template struct ThermoOsmosisModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
