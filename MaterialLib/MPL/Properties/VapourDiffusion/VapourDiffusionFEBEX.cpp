/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 5, 2021, 3:49 PM
 */

#include "VapourDiffusionFEBEX.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
PropertyDataType VapourDiffusionFEBEX::value(
    const VariableArray& variable_array,
    const ParameterLib::SpatialPosition& /*pos*/, const double /*t*/,
    const double /*dt*/) const
{
    const double T = variable_array.temperature;

    return base_diffusion_coefficient_ *
           std::pow(T / MaterialLib::PhysicalConstant::CelsiusZeroInKelvin,
                    exponent_);
}

PropertyDataType VapourDiffusionFEBEX::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double T = variable_array.temperature;

    if (variable == Variable::temperature)
    {
        return exponent_ * base_diffusion_coefficient_ *
               std::pow(T / MaterialLib::PhysicalConstant::CelsiusZeroInKelvin,
                        exponent_ - 1.0) /
               MaterialLib::PhysicalConstant::CelsiusZeroInKelvin;
    }
    if (variable == Variable::liquid_saturation)
    {
        return 0.0;
    }

    OGS_FATAL(
        "VapourDiffusionFEBEX::dValue is implemented for derivatives with "
        "respect to temperature or saturation only.");
}

}  // namespace MaterialPropertyLib
