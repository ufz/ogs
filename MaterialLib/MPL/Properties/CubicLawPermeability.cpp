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
CubicLawPermeability::CubicLawPermeability(
    std::string name, ParameterLib::Parameter<double> const* b)
    : _b(b)
{
    name_ = std::move(name);
}

PropertyDataType CubicLawPermeability::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    double const aperture_m =
        _b ? (*_b)(t, pos)[0] : variable_array.fracture_aperture;

    return aperture_m * aperture_m / 12;
}

PropertyDataType CubicLawPermeability::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (variable != Variable::fracture_aperture || _b)
    {
        return 0.0;
    }

    return variable_array.fracture_aperture / 6.0;
}
}  // namespace MaterialPropertyLib
