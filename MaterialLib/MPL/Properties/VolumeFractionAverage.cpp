/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 20, 2022
 */

#include "VolumeFractionAverage.h"

#include <iostream>

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
    for (auto& phase : phases)
    {
        if (phase == nullptr)
        {
            OGS_FATAL(
                "One of the required phases (AqueousLiquid/FrozenLiquid/Solid) "
                "does not exist!");
        }
        std::string phase_name = phase->name;

        if (!phase->hasProperty(prop_type_))
        {
            OGS_FATAL(
                "The phase '{}' does not have the required property '{}'!",
                phase_name, property_enum_to_string[prop_type_]);
        }
        auto const& property = phase->property(prop_type_);
        if (phase_name == "AqueousLiquid")
        {
            properties_.liquid = &property;
        }
        else if (phase_name == "FrozenLiquid")
        {
            properties_.frozen = &property;
        }
        else if (phase_name == "Solid")
        {
            properties_.porous = &property;
        }
    }
}

PropertyDataType VolumeFractionAverage::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    auto const& medium = *std::get<Medium*>(scale_);
    auto const& porosity = medium[PropertyType::porosity];
    auto const& fraction = medium[PropertyType::volume_fraction];

    // get frozen pore volume fraction, and porosity
    auto const pfr =
        std::get<double>(fraction.value(variable_array, pos, t, dt));
    auto const phi =
        std::get<double>(porosity.value(variable_array, pos, t, dt));

    auto const prop_value_liquid =
        std::get<double>(properties_.liquid->value(variable_array, pos, t, dt));
    auto const prop_value_frozen =
        std::get<double>(properties_.frozen->value(variable_array, pos, t, dt));
    auto const prop_value_porous =
        std::get<double>(properties_.porous->value(variable_array, pos, t, dt));

    auto average_value = (phi - pfr) * prop_value_liquid +
                         pfr * prop_value_frozen +
                         (1 - phi) * prop_value_porous;

    return average_value;
}

PropertyDataType VolumeFractionAverage::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    (void)variable;
    assert((variable == Variable::temperature) &&
           "VolumeFractionAverage::dvalue is implemented for "
           "derivatives with respect to temperature only.");

    auto const& medium = *std::get<Medium*>(scale_);
    auto const& fraction = medium[PropertyType::volume_fraction];

    auto const dpfr_dT = std::get<double>(
        fraction.dValue(variable_array, Variable::temperature, pos, t, dt));

    double prop_value_liquid =
        std::get<double>(properties_.liquid->value(variable_array, pos, t, dt));
    double prop_value_frozen =
        std::get<double>(properties_.frozen->value(variable_array, pos, t, dt));

    return (prop_value_frozen - prop_value_liquid) * dpfr_dT;
}
}  // namespace MaterialPropertyLib
