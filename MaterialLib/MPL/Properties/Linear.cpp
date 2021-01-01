/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <numeric>

#include "MaterialLib/MPL/Properties/Linear.h"

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
                               ParameterLib::SpatialPosition const& /*pos*/,
                               double const /*t*/, double const /*dt*/) const
{
    auto calculate_linearized_ratio = [&variable_array](
                                          double const initial_linearized_ratio,
                                          auto const& iv) {
        return initial_linearized_ratio +
               std::get<double>(iv.slope) *
                   (std::get<double>(
                        variable_array[static_cast<int>(iv.type)]) -
                    std::get<double>(iv.reference_condition));
    };

    double const linearized_ratio_to_reference_value =
        std::accumulate(independent_variables_.begin(),
                        independent_variables_.end(),
                        1.0,
                        calculate_linearized_ratio);

    return std::get<double>(value_) * linearized_ratio_to_reference_value;
}

PropertyDataType Linear::dValue(VariableArray const& /*variable_array*/,
                                Variable const primary_variable,
                                ParameterLib::SpatialPosition const& /*pos*/,
                                double const /*t*/, double const /*dt*/) const
{
    auto const independent_variable =
        std::find_if(independent_variables_.begin(),
                     independent_variables_.end(),
                     [&primary_variable](auto const& iv) -> bool {
                         return iv.type == primary_variable;
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
