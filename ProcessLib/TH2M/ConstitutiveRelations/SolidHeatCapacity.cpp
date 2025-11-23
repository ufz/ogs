/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SolidHeatCapacity.h"

namespace ProcessLib::TH2M
{
namespace ConstitutiveRelations
{

void SolidHeatCapacityModel::eval(
    SpaceTimeData const& x_t,
    MediaData const& media_data,
    TemperatureData const& T_data,
    SolidHeatCapacityData& solid_heat_capacity) const
{
    namespace MPL = MaterialPropertyLib;

    MPL::VariableArray variables;
    variables.temperature = T_data.T;

    auto const& mpl_cpS = media_data.specific_heat_capacity_solid;

    *solid_heat_capacity =
        mpl_cpS.template value<double>(variables, x_t.x, x_t.t, x_t.dt);
}

}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::TH2M
