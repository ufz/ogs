// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "VolumeFractionAverage.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/PropertyType.h"

namespace MaterialPropertyLib
{
VolumeFractionAverage::VolumeFractionAverage(std::string name)
{
    // get the corresponding property's name
    name_ = std::move(name);

    prop_type_ = convertStringToProperty(name_);
}

void VolumeFractionAverage::checkScale() const
{
    if (!std::holds_alternative<Medium*>(scale_))
    {
        OGS_FATAL(
            "The property 'VolumeFractionAverage' is "
            "implemented on the 'medium' scale only.");
    }
}

void VolumeFractionAverage::setProperties(
    std::vector<std::unique_ptr<Phase>> const& phases)
{
    // run over phases, identify them and get properties
    for (auto const& phase : phases)
    {
        if (phase == nullptr)
        {
            OGS_FATAL(
                "One of the required phases (AqueousLiquid/FrozenLiquid/Solid) "
                "does not exist!");
        }

        if (!phase->hasProperty(prop_type_))
        {
            OGS_FATAL(
                "The phase '{}' does not have the required property '{}'!",
                toString(phase->phaseName),
                property_enum_to_string[prop_type_]);
        }
        auto const& property = phase->property(prop_type_);
        if (phase->phaseName == PhaseName::AqueousLiquid)
        {
            properties_.liquid = &property;
        }
        else if (phase->phaseName == PhaseName::FrozenLiquid)
        {
            properties_.frozen = &property;
        }
        else if (phase->phaseName == PhaseName::Solid)
        {
            properties_.porous = &property;
        }
    }

    // Consistency check performed once at setup (scale_ is already set to the
    // medium by setScale() before setProperties() runs), instead of on every
    // value()/dValue() call: a FrozenLiquid phase providing this property and
    // the medium's frozen_liquid_saturation property must either both be
    // present or both be absent.
    auto const& medium = *std::get<Medium*>(scale_);
    if ((properties_.frozen != nullptr) !=
        medium.hasProperty(PropertyType::frozen_liquid_saturation))
    {
        OGS_FATAL(
            "In the 'VolumeFractionAverage' property '{}', a FrozenLiquid "
            "phase "
            "and the medium's 'frozen_liquid_saturation' property must be "
            "provided together. Found FrozenLiquid phase property: {}, medium "
            "frozen_liquid_saturation property: {}.",
            name_, properties_.frozen != nullptr,
            medium.hasProperty(PropertyType::frozen_liquid_saturation));
    }
}

PropertyDataType VolumeFractionAverage::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    auto const& medium = *std::get<Medium*>(scale_);
    auto const& porosity = medium[PropertyType::porosity];

    auto const phi =
        std::get<double>(porosity.value(variable_array, pos, t, dt));

    double S_fr = 0;
    double prop_value_frozen = 0;

    // get frozen liquid (ice) saturation, i.e. the fraction of the pore space
    // occupied by ice (the frozen-phase/saturation consistency is enforced once
    // in setProperties())
    if (medium.hasProperty(PropertyType::frozen_liquid_saturation))
    {
        auto const& saturation = medium[PropertyType::frozen_liquid_saturation];
        S_fr = std::get<double>(saturation.value(variable_array, pos, t, dt));
        prop_value_frozen = std::get<double>(
            properties_.frozen->value(variable_array, pos, t, dt));
    }
    auto const prop_value_liquid =
        std::get<double>(properties_.liquid->value(variable_array, pos, t, dt));
    auto const prop_value_porous =
        std::get<double>(properties_.porous->value(variable_array, pos, t, dt));

    return phi * (1.0 - S_fr) * prop_value_liquid +
           phi * S_fr * prop_value_frozen + (1 - phi) * prop_value_porous;
}

PropertyDataType VolumeFractionAverage::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    if (variable != Variable::temperature)
    {
        OGS_FATAL(
            "VolumeFractionAverage::dValue is implemented for derivatives with "
            "respect to temperature only.");
    }

    auto const& medium = *std::get<Medium*>(scale_);
    auto const& porosity = medium[PropertyType::porosity];
    auto const phi =
        std::get<double>(porosity.value(variable_array, pos, t, dt));

    double dS_fr_dT = 0;
    double prop_value_frozen = 0;

    if (medium.hasProperty(PropertyType::frozen_liquid_saturation))
    {
        auto const& saturation = medium[PropertyType::frozen_liquid_saturation];
        dS_fr_dT = std::get<double>(saturation.dValue(
            variable_array, Variable::temperature, pos, t, dt));
        prop_value_frozen = std::get<double>(
            properties_.frozen->value(variable_array, pos, t, dt));
    }

    double prop_value_liquid =
        std::get<double>(properties_.liquid->value(variable_array, pos, t, dt));

    // Only the frozen liquid saturation S_fr is treated as temperature
    // dependent here; the porosity and the phase properties are assumed to be
    // temperature independent, so their temperature derivatives are omitted.
    return phi * (prop_value_frozen - prop_value_liquid) * dS_fr_dT;
}
}  // namespace MaterialPropertyLib
