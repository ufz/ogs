// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

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
    double const pressure = variable_array.liquid_phase_pressure;
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
