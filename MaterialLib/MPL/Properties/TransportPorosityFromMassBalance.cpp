/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TransportPorosityFromMassBalance.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
void TransportPorosityFromMassBalance::checkScale() const
{
    if (!std::holds_alternative<Phase*>(scale_))
    {
        OGS_FATAL(
            "The property 'TransportPorosityFromMassBalance' is "
            "implemented on the 'phase' scales only.");
    }
    auto const phase = std::get<Phase*>(scale_);
    if (phase->name != "Solid")
    {
        OGS_FATAL(
            "The property 'TransportPorosityFromMassBalance' must be given in "
            "the 'Solid' phase, not in '{:s}' phase.",
            phase->name);
    }
}

PropertyDataType TransportPorosityFromMassBalance::value(
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

    double const e = std::get<double>(
        variable_array[static_cast<int>(Variable::volumetric_strain)]);
    double const e_prev = std::get<double>(
        variable_array_prev[static_cast<int>(Variable::volumetric_strain)]);
    double const delta_e = e - e_prev;

    double const p_eff = std::get<double>(
        variable_array[static_cast<int>(Variable::effective_pore_pressure)]);
    double const p_eff_prev =
        std::get<double>(variable_array_prev[static_cast<int>(
            Variable::effective_pore_pressure)]);
    double const delta_p_eff = p_eff - p_eff_prev;

    double const phi =
        std::get<double>(variable_array[static_cast<int>(Variable::porosity)]);

    double const phi_tr_prev = std::get<double>(
        variable_array_prev[static_cast<int>(Variable::transport_porosity)]);

    double const w = delta_e + delta_p_eff * beta_SR;
    return std::clamp(phi_tr_prev + (alpha_b - phi) * w, phi_min_, phi_max_);
}

PropertyDataType TransportPorosityFromMassBalance::value(
    VariableArray const& /*variable_array*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL(
        "TransportPorosityFromMassBalance value call requires previous time "
        "step values.");
}

PropertyDataType TransportPorosityFromMassBalance::dValue(
    VariableArray const& /*variable_array*/,
    Variable const /*primary_variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL(
        "TransportPorosityFromMassBalance derivatives are not implemented.");
}

}  // namespace MaterialPropertyLib
