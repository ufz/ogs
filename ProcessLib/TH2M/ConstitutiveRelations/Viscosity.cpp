/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Viscosity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{
void ViscosityModel::eval(SpaceTimeData const& x_t, MediaData const& media_data,
                          TemperatureData const& T_data,
                          MassMoleFractionsData const& mass_mole_fractions_data,
                          ViscosityData& viscosity_data) const
{
    MaterialPropertyLib::VariableArray variables;

    variables.temperature = T_data.T;
    variables.molar_fraction = mass_mole_fractions_data.xnCG;

    auto const& liquid_phase = media_data.liquid;
    auto const& gas_phase = media_data.gas;

    viscosity_data.mu_GR =
        gas_phase[MaterialPropertyLib::PropertyType::viscosity]
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);

    viscosity_data.mu_LR =
        liquid_phase[MaterialPropertyLib::PropertyType::viscosity]
            .template value<double>(variables, x_t.x, x_t.t, x_t.dt);
}
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
