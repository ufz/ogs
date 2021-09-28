/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "BishopsSaturationCutoff.h"

namespace MaterialPropertyLib
{
BishopsSaturationCutoff::BishopsSaturationCutoff(std::string name,
                                                 double const cutoff_value)
    : S_L_max_(cutoff_value)
{
    name_ = std::move(name);
}

void BishopsSaturationCutoff::checkScale() const
{
    if (!std::holds_alternative<Medium*>(scale_))
    {
        OGS_FATAL(
            "The property 'BishopsSaturationCutoff' is implemented on the "
            "'media' scale only.");
    }
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
    if (variable != Variable::liquid_saturation)
    {
        OGS_FATAL(
            "BishopsSaturationCutoff::dValue is implemented for derivatives "
            "with respect to liquid saturation only.");
    }

    return 0.;
}
}  // namespace MaterialPropertyLib
