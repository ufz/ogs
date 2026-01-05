// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MaterialLib/MPL/Properties/Exponential.h"

#include <cmath>

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

PropertyDataType Exponential::value(VariableArray const& variable_array,
                                    ParameterLib::SpatialPosition const& pos,
                                    double const t, double const /*dt*/) const
{
    double v = 0.0;
    if (Variable const* const variable =
            std::get_if<Variable>(&exponent_data_.type))
    {
        v = std::get<double>(variable_array[*variable]);
    }
    else if (auto* str_ptr = std::get_if<std::string>(&exponent_data_.type))
    {
        if (*str_ptr == "t")
            v = t;
        else if (*str_ptr == "x")
            v = pos.getCoordinates().value()[0];
        else if (*str_ptr == "y")
            v = pos.getCoordinates().value()[1];
        else if (*str_ptr == "z")
            v = pos.getCoordinates().value()[2];
        else
            OGS_FATAL(
                "Unknown independent_variable {:s} for exponential property.",
                *str_ptr)
    }
    else
    {
        OGS_FATAL(
            "Could not convert independent_variable neither to a Variable nor "
            "to a std::string.");
    }
    auto const f = std::get<double>(exponent_data_.factor);

    return offset_ + std::get<double>(value_) * std::exp(f * v);
}

PropertyDataType Exponential::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    Variable const* const independent_variable =
        std::get_if<Variable>(&exponent_data_.type);
    if (independent_variable == nullptr)
    {
        return 0.;
    }
    if (*independent_variable != variable)
    {
        return 0.;
    }

    auto const f = std::get<double>(exponent_data_.factor);
    auto const v = std::get<double>(variable_array[*independent_variable]);

    return std::get<double>(value_) * f * std::exp(f * v);
}

PropertyDataType Exponential::d2Value(
    VariableArray const& variable_array, Variable const pv1, Variable const pv2,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    Variable const* const independent_variable =
        std::get_if<Variable>(&exponent_data_.type);
    if (independent_variable == nullptr)
    {
        return 0.;
    }
    if (*independent_variable != pv1 && *independent_variable != pv2)
    {
        return 0.;
    }

    auto const f = std::get<double>(exponent_data_.factor);
    auto const v = std::get<double>(variable_array[*independent_variable]);

    return std::get<double>(value_) * f * f * std::exp(f * v);
}

}  // namespace MaterialPropertyLib
