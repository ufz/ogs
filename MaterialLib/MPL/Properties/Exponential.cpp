/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
                         PropertyDataType const& property_reference_value,
                         ExponentData const& v)
    : exponent_data_(v)
{
    name_ = std::move(name);
    value_ = property_reference_value;
}

PropertyDataType Exponential::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    return std::get<double>(value_) *
           std::exp(
               -std::get<double>(exponent_data_.factor) *
               (std::get<double>(
                    variable_array[static_cast<int>(exponent_data_.type)]) -
                std::get<double>(exponent_data_.reference_condition)));
}

PropertyDataType Exponential::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    return exponent_data_.type == primary_variable
               ? -std::get<double>(value_) *
                     std::get<double>(exponent_data_.factor) *
                     std::exp(
                         -std::get<double>(exponent_data_.factor) *
                         (std::get<double>(variable_array[static_cast<int>(
                              exponent_data_.type)]) -
                          std::get<double>(exponent_data_.reference_condition)))
               : decltype(value_){};
}

PropertyDataType Exponential::d2Value(
    VariableArray const& variable_array, Variable const pv1, Variable const pv2,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    return exponent_data_.type == pv1 && exponent_data_.type == pv2
               ? std::get<double>(value_) *
                     boost::math::pow<2>(
                         std::get<double>(exponent_data_.factor)) *
                     std::exp(
                         -std::get<double>(exponent_data_.factor) *
                         (std::get<double>(variable_array[static_cast<int>(
                              exponent_data_.type)]) -
                          std::get<double>(exponent_data_.reference_condition)))
               : decltype(value_){};
}

}  // namespace MaterialPropertyLib
