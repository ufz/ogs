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
    ParameterLib::Parameter<double> const& parameter)
    : _parameter(parameter)
{
}

PropertyDataType ParameterProperty::value(
    VariableArray const& /*variable_array*/,
    ParameterLib::SpatialPosition const& pos,
    double const t) const
{
    auto const& values = _parameter(t, pos);
    switch (values.size())
    {
        case 1:
        {
            // scalar
            return values[0];
        }
        case 2:
        {
            // Pair
            return Pair{values[0], values[1]};
        }
        case 3:
        {
            // Vector
            return Vector{values[0], values[1], values[2]};
        }
        case 4:
        {
            // Tensor
            return Tensor2d{values[0], values[1], values[2], values[3]};
        }
        case 6:
        {
            // Symmetric Tensor
            return SymmTensor{values[0], values[1], values[2],
                              values[3], values[4], values[5]};
        }
        case 9:
        {
            // Tensor
            return Tensor{values[0], values[1], values[2], values[3], values[4],
                          values[5], values[6], values[7], values[8]};
        }

        default:
        {
            OGS_FATAL(
                "Creation of a parameterized property with %i components is "
                "not implemented.",
                values.size());
        }
    }
}

PropertyDataType ParameterProperty::dValue(
    VariableArray const& /*variable_array*/,
    Variable const /*primary_variable*/,
    ParameterLib::SpatialPosition const& /*pos*/,
    double const /*t*/) const
{
    return double{};
}

PropertyDataType ParameterProperty::d2Value(
    VariableArray const& /*variable_array*/,
    Variable const /*pv1*/,
    Variable const /*pv2*/,
    ParameterLib::SpatialPosition const& /*pos*/,
    double const /*t*/) const
{
    return double{};
}

}  // namespace MaterialPropertyLib
