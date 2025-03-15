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
        double x = 0.0;
        if (auto* var_ptr = std::get_if<Variable>(&iv.type))
        {
            x = std::get<double>(variable_array[*var_ptr]);
        }
        else if (auto* str_ptr = std::get_if<std::string>(&iv.type))
        {
            if (*str_ptr == "t")
                x = t;
            else if (*str_ptr == "x")
                x = pos.getCoordinates().value()[0];
            else if (*str_ptr == "y")
                x = pos.getCoordinates().value()[1];
            else if (*str_ptr == "z")
                x = pos.getCoordinates().value()[2];
            else
                OGS_FATAL(
                    "Unknown independent_variable {:s} for curve property.",
                    *str_ptr)
        }
        else
        {
            OGS_FATAL(
                "Could not convert independent_variable neither to a Variable "
                "nor to a std::string.");
        }

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

PropertyDataType Linear::dValue(VariableArray const& /*variable_array*/,
                                Variable const variable,
                                ParameterLib::SpatialPosition const& /*pos*/,
                                double const /*t*/, double const /*dt*/) const
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

    return independent_variable != independent_variables_.end()
               ? std::get<double>(value_) *
                     std::get<double>(independent_variable->slope)
               : decltype(value_){};
}

PropertyDataType Linear::d2Value(VariableArray const& /*variable_array*/,
                                 Variable const /*pv1*/, Variable const /*pv2*/,
                                 ParameterLib::SpatialPosition const& /*pos*/,
                                 double const /*t*/, double const /*dt*/) const
{
    return decltype(value_){};
}

}  // namespace MaterialPropertyLib
