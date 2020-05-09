/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "MaterialLib/MPL/Properties/ParameterProperty.h"

namespace MaterialPropertyLib
{
ParameterProperty::ParameterProperty(
    std::string name, ParameterLib::Parameter<double> const& parameter)
    : parameter_(parameter)
{
    name_ = std::move(name);
}

PropertyDataType ParameterProperty::value(
    VariableArray const& /*variable_array*/,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    return fromVector(parameter_(t, pos));
}

PropertyDataType ParameterProperty::dValue(
    VariableArray const& /*variable_array*/,
    Variable const /*primary_variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    return double{};
}

PropertyDataType ParameterProperty::d2Value(
    VariableArray const& /*variable_array*/, Variable const /*pv1*/,
    Variable const /*pv2*/, ParameterLib::SpatialPosition const& /*pos*/,
    double const /*t*/, double const /*dt*/) const
{
    return double{};
}

}  // namespace MaterialPropertyLib
