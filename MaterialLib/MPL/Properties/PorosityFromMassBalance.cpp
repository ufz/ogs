/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "PorosityFromMassBalance.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
void PorosityFromMassBalance::checkScale() const
{
    if (!std::holds_alternative<Phase*>(scale_))
    {
        OGS_FATAL(
            "The property 'PorosityFromMassBalance' is "
            "implemented on the 'phase' scales only.");
    }
    auto const phase = std::get<Phase*>(scale_);
    if (phase->name != "Solid")
    {
        OGS_FATAL(
            "The property 'PorosityFromMassBalance' must be given in the "
            "'Solid' phase, not in '{:s}' phase.",
            phase->name);
    }
}

PropertyDataType PorosityFromMassBalance::value(
    VariableArray const& /*variable_array*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL(
        "PorosityFromMassBalance value call requires previous time step "
        "values.");
}

PropertyDataType PorosityFromMassBalance::value(
    VariableArray const& variable_array,
    VariableArray const& variable_array_prev,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    double const beta_SR = std::get<double>(
        variable_array[static_cast<int>(Variable::grain_compressibility)]);
    auto const alpha_b =
        std::get<Phase*>(scale_)
            ->property(PropertyType::biot_coefficient)
            .template value<double>(variable_array, pos, t, dt);

    double const e_dot = std::get<double>(
        variable_array[static_cast<int>(Variable::volumetric_strain_rate)]);

    double const p_eff_dot = std::get<double>(variable_array[static_cast<int>(
        Variable::effective_pore_pressure_rate)]);

    double const phi_prev = std::get<double>(
        variable_array_prev[static_cast<int>(Variable::porosity)]);

    double const w = dt * (e_dot + p_eff_dot * beta_SR);
    return std::clamp((phi_prev + alpha_b * w) / (1 + w), phi_min_, phi_max_);
}

PropertyDataType PorosityFromMassBalance::dValue(
    VariableArray const& /*variable_array*/,
    Variable const /*primary_variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL("PorosityFromMassBalance derivatives are not implemented.");
}

}  // namespace MaterialPropertyLib
