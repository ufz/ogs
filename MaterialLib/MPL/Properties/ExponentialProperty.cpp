/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include <boost/math/special_functions/pow.hpp>
#include <cmath>

#include "MaterialLib/MPL/Properties/ExponentialProperty.h"

namespace MaterialPropertyLib
{
ExponentialProperty::ExponentialProperty(
    PropertyDataType const& property_reference_value, ExponentData const& v)
    : _exponent_data(v)
{
    _value = property_reference_value;
}

PropertyDataType ExponentialProperty::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/,
    double const /*t*/) const
{
    return std::get<double>(_value) *
           std::exp(
               -std::get<double>(_exponent_data.factor) *
               (std::get<double>(
                    variable_array[static_cast<int>(_exponent_data.type)]) -
                std::get<double>(_exponent_data.reference_condition)));
}

PropertyDataType ExponentialProperty::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/) const
{
    return _exponent_data.type == primary_variable
               ? -std::get<double>(_value) *
                     std::get<double>(_exponent_data.factor) *
                     std::exp(
                         -std::get<double>(_exponent_data.factor) *
                         (std::get<double>(variable_array[static_cast<int>(
                              _exponent_data.type)]) -
                          std::get<double>(_exponent_data.reference_condition)))
               : decltype(_value){};
}

PropertyDataType ExponentialProperty::d2Value(
    VariableArray const& variable_array,
    Variable const pv1,
    Variable const pv2) const
{
    return _exponent_data.type == pv1 && _exponent_data.type == pv2
               ? std::get<double>(_value) *
                     boost::math::pow<2>(
                         std::get<double>(_exponent_data.factor)) *
                     std::exp(
                         -std::get<double>(_exponent_data.factor) *
                         (std::get<double>(variable_array[static_cast<int>(
                              _exponent_data.type)]) -
                          std::get<double>(_exponent_data.reference_condition)))
               : decltype(_value){};
}

}  // namespace MaterialPropertyLib
