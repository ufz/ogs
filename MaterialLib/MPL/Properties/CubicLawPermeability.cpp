/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CubicLawPermeability.h"

namespace MaterialPropertyLib
{
PropertyDataType CubicLawPermeability::value(
    VariableArray const& /*variable_array*/,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    double const aperture_m = _b(t, pos)[0];
    return aperture_m * aperture_m / 12;
}

PropertyDataType CubicLawPermeability::dValue(
    VariableArray const& /*variable_array*/, Variable const /*variable*/,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    OGS_FATAL("CubicLawPermeability::dValue is not implemented.");
}
}  // namespace MaterialPropertyLib
