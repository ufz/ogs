/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BishopsSaturationCutoff.h"

namespace MaterialPropertyLib
{
BishopsSaturationCutoff::BishopsSaturationCutoff(double const cutoff_value)
    : S_L_max_(cutoff_value)
{
}

void BishopsSaturationCutoff::setScale(
    std::variant<Medium*, Phase*, Component*> scale_pointer)
{
    if (!std::holds_alternative<Medium*>(scale_pointer))
    {
        OGS_FATAL(
            "The property 'BishopsSaturationCutoff' is implemented on the "
            "'media' scale only.");
    }
    medium_ = std::get<Medium*>(scale_pointer);
}

PropertyDataType BishopsSaturationCutoff::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    auto const S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    return S_L < S_L_max_ ? 0. : 1.;
}

PropertyDataType BishopsSaturationCutoff::dValue(
    VariableArray const& /*variable_array*/, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)variable;
    assert(
        (variable == Variable::liquid_saturation) &&
        "BishopsSaturationCutoff::dvalue is implemented for derivatives with "
        "respect to liquid saturation only.");

    return 0.;
}
}  // namespace MaterialPropertyLib
