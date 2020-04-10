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
void PorosityFromMassBalance::setScale(
    std::variant<Medium*, Phase*, Component*> scale_pointer)
{
    if (std::holds_alternative<Phase*>(scale_pointer))
    {
        _phase = std::get<Phase*>(scale_pointer);
        if (_phase->name != "Solid")
        {
            OGS_FATAL(
                "The property 'PorosityFromMassBalance' must be "
                "given in the 'Solid' phase, not in '{:s}' phase.",
                _phase->name.c_str());
        }
    }
    else
    {
        OGS_FATAL(
            "The property 'PorosityFromMassBalance' is "
            "implemented on the 'phase' scales only.");
    }
}

PropertyDataType PorosityFromMassBalance::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos,
    double const t, double const dt) const
{
    double const K_SR = _phase->property(PropertyType::bulk_modulus)
                            .template value<double>(variable_array, pos, t, dt);
    auto const alpha_b =
        _phase->property(PropertyType::biot_coefficient)
            .template value<double>(variable_array, pos, t, dt);

    double const e_dot = std::get<double>(
        variable_array[static_cast<int>(Variable::volumetric_strain_rate)]);

    double const p_eff_dot = std::get<double>(variable_array[static_cast<int>(
        Variable::effective_pore_pressure_rate)]);

    double const phi = std::get<double>(
        variable_array[static_cast<int>(Variable::porosity)]);

    double const w = dt * (e_dot + p_eff_dot / K_SR);
    return std::clamp((phi + alpha_b * w) / (1 + w), _phi_min, _phi_max);
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
