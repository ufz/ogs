/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DarcyLaw.h"

#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void DarcyLawModel<DisplacementDim>::eval(
    CapillaryPressureData<DisplacementDim> const& p_cap_data,
    LiquidDensityData const& rho_L_data,
    LiquidViscosityData const& mu_L_data,
    PermeabilityData<DisplacementDim> const& perm_data,
    ThermoOsmosisData<DisplacementDim> const& th_osmosis_data,
    DarcyLawData<DisplacementDim>& out) const
{
    out.v_darcy = perm_data.Ki / mu_L_data.viscosity *
                      (perm_data.k_rel *
                       (p_cap_data.grad_p_cap + rho_L_data.rho_LR * b_)) +
                  th_osmosis_data.seepage_velocity_contribution;
}

template struct DarcyLawModel<2>;
template struct DarcyLawModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
