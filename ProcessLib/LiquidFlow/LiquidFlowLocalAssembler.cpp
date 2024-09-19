/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on September 30, 2024, 11:54 AM
 */

#include "LiquidFlowLocalAssembler.h"

namespace ProcessLib
{
namespace LiquidFlow
{
std::tuple<double, double> getFluidDensityAndViscosity(
    double const t, double const dt, ParameterLib::SpatialPosition const& pos,
    MaterialPropertyLib::Phase const& fluid_phase,
    MaterialPropertyLib::VariableArray& vars)
{
    // Compute density:
    // Quick workaround: If fluid density is described as ideal gas, then
    // the molar mass must be passed to the MPL::IdealGasLaw via the
    // variable_array and the fluid must have the property
    // MPL::PropertyType::molar_mass. For other density models (e.g.
    // Constant), it is not mandatory to specify the molar mass.
    if (fluid_phase.hasProperty(MaterialPropertyLib::PropertyType::molar_mass))
    {
        vars.molar_mass =
            fluid_phase.property(MaterialPropertyLib::PropertyType::molar_mass)
                .template value<double>(vars, pos, t, dt);
    }

    auto const fluid_density =
        fluid_phase[MaterialPropertyLib::PropertyType::density]
            .template value<double>(vars, pos, t, dt);
    assert(fluid_density > 0.);
    vars.density = fluid_density;

    auto const viscosity =
        fluid_phase[MaterialPropertyLib::PropertyType::viscosity]
            .template value<double>(vars, pos, t, dt);

    return {fluid_density, viscosity};
}
}  // namespace LiquidFlow
}  // namespace ProcessLib