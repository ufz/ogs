/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MaterialLib/MPL/Properties/Linear.h"

#include <numeric>

namespace MaterialPropertyLib
{
double IndependentVariable::valueClamped(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos,
    double const t) const
{
    double x = valueUnclamped(variable_array, pos, t);

    if (min)
    {
        x = std::max(x, *min);
    }
    if (max)
    {
        x = std::min(x, *max);
    }

    return x;
}

double IndependentVariable::valueUnclamped(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos,
    double const t) const
{
    if (auto* var_ptr = std::get_if<Variable>(&type))
    {
        return std::get<double>(variable_array[*var_ptr]);
    }

    if (auto* str_ptr = std::get_if<std::string>(&type))
    {
        if (*str_ptr == "t")
            return t;

        if (*str_ptr == "x")
            return pos.getCoordinates().value()[0];

        if (*str_ptr == "y")
            return pos.getCoordinates().value()[1];

        if (*str_ptr == "z")
            return pos.getCoordinates().value()[2];

        OGS_FATAL("Unknown independent variable {:s} for a linear property.",
                  *str_ptr)
    }

    OGS_FATAL(
        "Could not convert independent_variable neither to a Variable nor "
        "to a std::string.");
}

Linear::Linear(std::string name,
               PropertyDataType const& property_reference_value,
               std::vector<IndependentVariable> const& vs)
    : independent_variables_(vs)
{
    name_ = std::move(name);
    value_ = property_reference_value;
}

PropertyDataType Linear::value(VariableArray const& variable_array,
                               ParameterLib::SpatialPosition const& pos,
                               double const t, double const /*dt*/) const
{
    auto calculate_linearized_ratio =
        [&variable_array, pos, t](double const initial_linearized_ratio,
                                  auto const& iv)
    {
        double const x = iv.valueClamped(variable_array, pos, t);

        return initial_linearized_ratio +
               std::get<double>(iv.slope) *
                   (x - std::get<double>(iv.reference_condition));
    };

    double const linearized_ratio_to_reference_value =
        std::accumulate(independent_variables_.begin(),
                        independent_variables_.end(),
                        1.0,
                        calculate_linearized_ratio);

    return std::get<double>(value_) * linearized_ratio_to_reference_value;
}

PropertyDataType Linear::dValue(VariableArray const& variable_array,
                                Variable const variable,
                                ParameterLib::SpatialPosition const& pos,
                                double const t, double const /*dt*/) const
{
    auto const independent_variable = std::find_if(
        independent_variables_.begin(),
        independent_variables_.end(),
        [&variable](auto const& iv) -> bool
        {
            if (auto const* var_ptr = std::get_if<Variable>(&iv.type))
            {
                return *var_ptr == variable;
            }
            return false;
        });

    auto const zero = decltype(value_){};
    if (independent_variable == independent_variables_.end())
    {
        return zero;
    }

    auto const& iv = *independent_variable;

    double const x = iv.valueUnclamped(variable_array, pos, t);

    if (iv.min && x < *iv.min)
    {
        return zero;
    }
    if (iv.max && x > *iv.max)
    {
        return zero;
    }

    return std::get<double>(value_) * std::get<double>(iv.slope);
}

PropertyDataType Linear::d2Value(VariableArray const& /*variable_array*/,
                                 Variable const /*pv1*/, Variable const /*pv2*/,
                                 ParameterLib::SpatialPosition const& /*pos*/,
                                 double const /*t*/, double const /*dt*/) const
{
    return decltype(value_){};
}

}  // namespace MaterialPropertyLib
