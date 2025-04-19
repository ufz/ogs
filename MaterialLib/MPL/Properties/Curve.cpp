/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Curve.h"

namespace MaterialPropertyLib
{
Curve::Curve(std::string name,
             StringOrVariable const independent_variable,
             MathLib::PiecewiseLinearInterpolation const& curve)
    : independent_variable_(independent_variable), curve_(curve)
{
    name_ = std::move(name);
}

PropertyDataType Curve::value(VariableArray const& variable_array,
                              ParameterLib::SpatialPosition const& pos,
                              double const t, double const /*dt*/) const
{
    double x = 0.0;
    if (Variable const* const independent_variable =
            std::get_if<Variable>(&independent_variable_))
    {
        x = std::get<double>(variable_array[*independent_variable]);
    }
    else if (std::string const* const str_ptr =
                 std::get_if<std::string>(&independent_variable_))
    {
        if (*str_ptr == "t")
        {
            x = t;
        }
        else if (*str_ptr == "x")
        {
            x = pos.getCoordinates().value()[0];
        }
        else if (*str_ptr == "y")
        {
            x = pos.getCoordinates().value()[1];
        }
        else if (*str_ptr == "z")
        {
            x = pos.getCoordinates().value()[2];
        }
        else
        {
            OGS_FATAL("Unknown independent_variable {:s} for curve property.",
                      *str_ptr)
        }
    }
    else
    {
        OGS_FATAL(
            "Could not convert independent_variable neither to a Variable nor "
            "to a std::string.");
    }
    return curve_.getValue(x);
}

PropertyDataType Curve::dValue(VariableArray const& variable_array,
                               Variable const variable,
                               ParameterLib::SpatialPosition const& /*pos*/,
                               double const /*t*/, double const /*dt*/) const
{
    Variable const* const independent_variable =
        std::get_if<Variable>(&independent_variable_);
    if (independent_variable == nullptr)
    {
        return 0.0;
    }
    if (variable != *independent_variable)
    {
        return 0.0;
    }
    auto const x = std::get<double>(variable_array[*independent_variable]);
    return curve_.getDerivative(x);
}
}  // namespace MaterialPropertyLib
