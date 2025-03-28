/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SolidDensity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{

void SolidDensityModel::eval(SpaceTimeData const& x_t,
                             MediaData const& media_data,
                             TemperatureData const& T_data,
                             SolidDensityData& solid_density_data) const
{
    MaterialPropertyLib::VariableArray variables;
    variables.temperature = T_data.T;

    auto const& mpl_solid_density =
        media_data.solid[MaterialPropertyLib::PropertyType::density];

    solid_density_data.rho_SR = mpl_solid_density.template value<double>(
        variables, x_t.x, x_t.t, x_t.dt);
}

void SolidDensityModel::dEval(
    SpaceTimeData const& x_t,
    MediaData const& media_data,
    TemperatureData const& T_data,
    SolidDensityDerivativeData& solid_density_d_data) const
{
    MaterialPropertyLib::VariableArray variables;
    variables.temperature = T_data.T;

    auto const& mpl_solid_density =
        media_data.solid[MaterialPropertyLib::PropertyType::density];

    solid_density_d_data.drho_SR_dT = mpl_solid_density.template dValue<double>(
        variables, MaterialPropertyLib::Variable::temperature, x_t.x, x_t.t,
        x_t.dt);
}

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
