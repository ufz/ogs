// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "HTFEM.h"

namespace ProcessLib
{
namespace HT
{
double evalEffectiveThermalExpansivity(
    double const t, double const dt, ParameterLib::SpatialPosition const& pos,
    MaterialPropertyLib::VariableArray const& vars,
    MaterialPropertyLib::Medium const& medium,
    MaterialPropertyLib::Phase const& liquid_phase,
    MaterialPropertyLib::Phase const& solid_phase,
    bool const has_solid_thermal_expansivity)
{
    double const dfluid_density_dT =
        liquid_phase.property(MaterialPropertyLib::PropertyType::density)
            .template dValue<double>(
                vars, MaterialPropertyLib::Variable::temperature, pos, t, dt);

    double const fluid_density = vars.density;
    double const porosity = vars.porosity;
    double const fluid_thermal_expansivity =
        -porosity * dfluid_density_dT / fluid_density;

    if (!has_solid_thermal_expansivity)
    {
        return fluid_thermal_expansivity;
    }

    double const linear_solid_thermal_expansivity =
        solid_phase
            .property(MaterialPropertyLib::PropertyType::thermal_expansivity)
            .template value<double>(vars, pos, t, dt);
    double const biot_coefficient =
        medium.property(MaterialPropertyLib::PropertyType::biot_coefficient)
            .template value<double>(vars, pos, t, dt);

    return fluid_thermal_expansivity + 3.0 * (biot_coefficient - porosity) *
                                           linear_solid_thermal_expansivity;
}
}  // namespace HT
}  // namespace ProcessLib
