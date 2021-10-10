/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MaterialLib/MPL/Properties/DupuitPermeability.h"

namespace MaterialPropertyLib
{
DupuitPermeability::DupuitPermeability(
    std::string name, ParameterLib::Parameter<double> const& parameter)
    : parameter_(parameter)
{
    name_ = std::move(name);
}

PropertyDataType DupuitPermeability::value(
    MaterialPropertyLib::VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    double const pressure = std::get<double>(variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::phase_pressure)]);
    auto const& permeability_values = parameter_(t, pos);

    auto const& permeability_variant = fromVector(permeability_values);
    PropertyDataType dupuit_permeability = std::visit(
        [&pressure](auto const& permeability_variant) -> PropertyDataType
        {
            using T = std::decay_t<decltype(permeability_variant)>;
            if constexpr (std::is_same_v<T, double>)
            {
                return pressure * permeability_variant;
            }
            else
            {
                return (pressure * permeability_variant).eval();
            }
        },
        permeability_variant);

    return dupuit_permeability;
}

}  // namespace MaterialPropertyLib
