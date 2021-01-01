/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <boost/math/special_functions/pow.hpp>
#include <cmath>

#include "MaterialLib/MPL/Properties/Exponential.h"

namespace MaterialPropertyLib
{
Exponential::Exponential(std::string name,
                         double const offset,
                         PropertyDataType const& property_reference_value,
                         ExponentData const& v)
    : exponent_data_(v), offset_(offset)
{
    name_ = std::move(name);
    auto const f = std::get<double>(exponent_data_.factor);
    auto const v0 = std::get<double>(exponent_data_.reference_condition);
    value_ = std::get<double>(property_reference_value) * std::exp(-f * v0);
}

PropertyDataType Exponential::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    auto const f = std::get<double>(exponent_data_.factor);
    auto const v =
        std::get<double>(variable_array[static_cast<int>(exponent_data_.type)]);

    return offset_ + std::get<double>(value_) * std::exp(f * v);
}

PropertyDataType Exponential::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (exponent_data_.type != primary_variable)
    {
        return 0.;
    }

    auto const f = std::get<double>(exponent_data_.factor);
    auto const v =
        std::get<double>(variable_array[static_cast<int>(exponent_data_.type)]);

    return std::get<double>(value_) * f * std::exp(f * v);
}

PropertyDataType Exponential::d2Value(
    VariableArray const& variable_array, Variable const pv1, Variable const pv2,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (exponent_data_.type != pv1 && exponent_data_.type != pv2)
    {
        return 0.;
    }

    auto const f = std::get<double>(exponent_data_.factor);
    auto const v =
        std::get<double>(variable_array[static_cast<int>(exponent_data_.type)]);

    return std::get<double>(value_) * f * f * std::exp(f * v);
}

}  // namespace MaterialPropertyLib
