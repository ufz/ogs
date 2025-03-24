/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"
#include "BaseLib/StrongType.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
using PureLiquidDensityData =
    BaseLib::StrongType<double, struct PureLiquidDensityTag>;

struct PureLiquidDensityModel
{
    void eval(SpaceTimeData const& x_t,
              MediaData const& media_data,
              GasPressureData const& p_GR,
              CapillaryPressureData const& p_cap,
              TemperatureData const& T_data,
              PureLiquidDensityData& out) const
    {
        MaterialPropertyLib::VariableArray variables;

        // primary variables
        auto const pGR = p_GR.pG;
        auto const pCap = p_cap.pCap;
        auto const T = T_data.T;
        variables.gas_phase_pressure = pGR;
        variables.capillary_pressure = pCap;
        variables.temperature = T;

        auto const& liquid_phase = media_data.liquid;

        // Water pressure is passed to the VariableArray in order to calculate
        // the water density.
        auto const pLR = pGR - pCap;
        variables.liquid_phase_pressure = pLR;

        // Concentration is initially zero to calculate the density of the pure
        // water phase, which is needed for the Kelvin-Laplace equation.
        variables.concentration = 0.;
        *out = liquid_phase.property(MaterialPropertyLib::PropertyType::density)
                   .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
    }
};

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
