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

#include "MaterialLib/MPL/Properties/LinearProperty.h"

namespace MaterialPropertyLib
{
LinearProperty::LinearProperty(PropertyDataType const& property_reference_value,
                               IndependentVariable const& v)
    : _independent_variable(v)
{
    _value = property_reference_value;
}

PropertyDataType LinearProperty::value(
    VariableArray const& variable_array) const
{
    return boost::get<double>(_value) +
           boost::get<double>(_independent_variable.slope) *
               (boost::get<double>(variable_array[_independent_variable.type]) -
                boost::get<double>(_independent_variable.reference_condition));
}

PropertyDataType LinearProperty::dValue(VariableArray const& /*variable_array*/,
                                        Variables const primary_variable) const
{
    return _independent_variable.type == primary_variable
               ? _independent_variable.slope
               : decltype(_value){};
}

PropertyDataType LinearProperty::d2Value(
    VariableArray const& /*variable_array*/,
    Variables const /*pv1*/,
    Variables const /*pv2*/) const
{
    return decltype(_value){};
}

}  // namespace MaterialPropertyLib
