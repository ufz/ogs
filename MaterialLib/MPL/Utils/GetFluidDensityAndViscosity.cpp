// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "GetFluidDensityAndViscosity.h"

#include "MaterialLib/MPL/Phase.h"
#include "MaterialLib/MPL/VariableType.h"
#include "ParameterLib/SpatialPosition.h"

namespace MaterialPropertyLib
{

double getFluidDensity(double const t, double const dt,
                       ParameterLib::SpatialPosition const& pos,
                       Phase const& fluid_phase, VariableArray& vars)
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

    return fluid_phase[MaterialPropertyLib::PropertyType::density]
        .template value<double>(vars, pos, t, dt);
}

std::tuple<double, double> getFluidDensityAndViscosity(
    double const t, double const dt, ParameterLib::SpatialPosition const& pos,
    MaterialPropertyLib::Phase const& fluid_phase,
    MaterialPropertyLib::VariableArray& vars)
{
    auto const fluid_density = getFluidDensity(t, dt, pos, fluid_phase, vars);

    assert(fluid_density > 0.);
    vars.density = fluid_density;

    auto const viscosity =
        fluid_phase[MaterialPropertyLib::PropertyType::viscosity]
            .template value<double>(vars, pos, t, dt);

    return {fluid_density, viscosity};
}

}  // namespace MaterialPropertyLib
