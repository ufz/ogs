// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "MaterialLib/MPL/Properties/Parameter.h"

namespace MaterialPropertyLib
{
Parameter::Parameter(std::string name,
                     ParameterLib::Parameter<double> const& parameter)
    : parameter_(parameter)
{
    name_ = std::move(name);
}

PropertyDataType Parameter::value(VariableArray const& /*variable_array*/,
                                  ParameterLib::SpatialPosition const& pos,
                                  double const t, double const /*dt*/) const
{
    return fromVector(parameter_(t, pos));
}

PropertyDataType Parameter::value(VariableArray const& /*variable_array*/,
                                  VariableArray const& /*variable_array_prev*/,
                                  ParameterLib::SpatialPosition const& pos,
                                  double const t, double const /*dt*/) const
{
    return fromVector(parameter_(t, pos));
}

PropertyDataType Parameter::dValue(VariableArray const& /*variable_array*/,
                                   Variable const /*variable*/,
                                   ParameterLib::SpatialPosition const& /*pos*/,
                                   double const /*t*/,
                                   double const /*dt*/) const
{
    return double{};
}

PropertyDataType Parameter::d2Value(
    VariableArray const& /*variable_array*/, Variable const /*pv1*/,
    Variable const /*pv2*/, ParameterLib::SpatialPosition const& /*pos*/,
    double const /*t*/, double const /*dt*/) const
{
    return double{};
}

}  // namespace MaterialPropertyLib
