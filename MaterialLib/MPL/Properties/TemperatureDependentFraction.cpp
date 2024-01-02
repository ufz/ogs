/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on May 20, 2022
 */

#include "TemperatureDependentFraction.h"

#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{
TemperatureDependentFraction::TemperatureDependentFraction(std::string name,
                                                           double const k,
                                                           double const T_c)
    : phase_change_shape_(k, T_c)
{
    name_ = std::move(name);
}

void TemperatureDependentFraction::checkScale() const
{
    if (!std::holds_alternative<Medium*>(scale_))
    {
        OGS_FATAL(
            "The property 'TemperatureDependantFraction' is "
            "implemented on the 'medium' scale only.");
    }
}

PropertyDataType TemperatureDependentFraction::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    auto const T = variable_array.temperature;

    auto const& medium = *std::get<Medium*>(scale_);
    auto const& porosity = medium[PropertyType::porosity];

    auto const phi =
        std::get<double>(porosity.value(variable_array, pos, t, dt));

    return phi * phase_change_shape_.value(T);
}

PropertyDataType TemperatureDependentFraction::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const dt) const
{
    (void)variable;
    assert((variable == Variable::temperature) &&
           "TemperatureDependantFraction::dvalue is implemented for "
           "derivatives with respect to temperature only.");

    auto const T = variable_array.temperature;

    auto const& medium = *std::get<Medium*>(scale_);
    auto const& porosity = medium[PropertyType::porosity];

    auto const phi =
        std::get<double>(porosity.value(variable_array, pos, t, dt));

    return phi * phase_change_shape_.dValue(T);
}

PropertyDataType TemperatureDependentFraction::d2Value(
    VariableArray const& variable_array, Variable const variable1,
    Variable const variable2, ParameterLib::SpatialPosition const& pos,
    double const t, double const dt) const
{
    (void)variable1;
    assert((variable1 == Variable::temperature) &&
           "TemperatureDependantFraction::d2value is implemented for "
           "derivatives with respect to temperature only.");

    (void)variable2;
    assert((variable2 == Variable::temperature) &&
           "TemperatureDependantFraction::d2value is implemented for "
           "derivatives with respect to temperature only.");

    auto const T = variable_array.temperature;

    auto const& medium = *std::get<Medium*>(scale_);
    auto const& porosity = medium[PropertyType::porosity];

    auto const phi =
        std::get<double>(porosity.value(variable_array, pos, t, dt));

    return phi * phase_change_shape_.d2Value(T);
}
}  // namespace MaterialPropertyLib
