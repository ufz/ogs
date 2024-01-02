/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
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
    if (!std::holds_alternative<Medium*>(scale_))
    {
        OGS_FATAL(
            "The property 'TransportPorosityFromMassBalance' is "
            "implemented on the 'medium' scales only.");
    }
}

PropertyDataType TransportPorosityFromMassBalance::value(
    VariableArray const& variable_array,
    VariableArray const& variable_array_prev,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    double const beta_SR = variable_array.grain_compressibility;
    auto const alpha_b =
        std::get<Medium*>(scale_)
            ->property(PropertyType::biot_coefficient)
            .template value<double>(variable_array, pos, t, dt);

    double const e = variable_array.volumetric_strain;
    double const e_prev = variable_array_prev.volumetric_strain;
    double const delta_e = e - e_prev;

    double const p_eff = variable_array.effective_pore_pressure;
    double const p_eff_prev = variable_array_prev.effective_pore_pressure;
    double const delta_p_eff = p_eff - p_eff_prev;

    double const phi = variable_array.porosity;

    double const phi_tr_prev = variable_array_prev.transport_porosity;

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
    VariableArray const& /*variable_array*/, Variable const /*variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL(
        "TransportPorosityFromMassBalance derivatives are not implemented.");
}

}  // namespace MaterialPropertyLib
