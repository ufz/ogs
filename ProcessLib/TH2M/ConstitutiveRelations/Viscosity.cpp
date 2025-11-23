/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
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

    viscosity_data.mu_GR = media_data.viscosity_gas.template value<double>(
        variables, x_t.x, x_t.t, x_t.dt);

    viscosity_data.mu_LR = media_data.viscosity_liquid.template value<double>(
        variables, x_t.x, x_t.t, x_t.dt);
}
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
