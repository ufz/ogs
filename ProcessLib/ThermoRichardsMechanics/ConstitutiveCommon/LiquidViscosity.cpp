/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "LiquidViscosity.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <int DisplacementDim>
void LiquidViscosityModel<DisplacementDim>::eval(
    SpaceTimeData const& x_t, MediaData const& media_data,
    LiquidDensityData const& rho_L_data,
    TemperatureData<DisplacementDim> const& T_data,
    LiquidViscosityData& out) const
{
    namespace MPL = MaterialPropertyLib;
    MPL::VariableArray variables;
    variables.temperature = T_data.T;
    variables.density = rho_L_data.rho_LR;

    out.viscosity =
        media_data.liquid.property(MPL::PropertyType::viscosity)
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
}

template struct LiquidViscosityModel<2>;
template struct LiquidViscosityModel<3>;
}  // namespace ProcessLib::ThermoRichardsMechanics
