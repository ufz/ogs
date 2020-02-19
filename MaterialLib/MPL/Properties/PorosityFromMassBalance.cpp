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

    double const p_dot = std::get<double>(
        variable_array[static_cast<int>(Variable::phase_pressure_rate)]);

    double const phi = std::get<double>(
        variable_array[static_cast<int>(Variable::porosity)]);

    double const w = dt * (e_dot + p_dot / K_SR);
    return (phi + alpha_b * w) / (1 + w);
}

PropertyDataType PorosityFromMassBalance::dValue(
    VariableArray const& /*variable_array*/, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::liquid_saturation) &&
           "PorosityFromMassBalance::dValue is implemented for "
           "derivatives with respect to liquid saturation only.");
    return 0;
}

}  // namespace MaterialPropertyLib
