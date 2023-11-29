/**
 * \file
 * \copyright
 * Copyright (c) 2012-2023, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 7, 2021, 9:14 AM
 */

#include "VapourDiffusionPMQ.h"

#include <algorithm>
#include <cmath>

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/VariableType.h"
#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
PropertyDataType VapourDiffusionPMQ::value(
    const VariableArray& variable_array,
    const ParameterLib::SpatialPosition& /*pos*/, const double /*t*/,
    const double /*dt*/) const
{
    const double S_L = std::clamp(variable_array.liquid_saturation, 0.0, 1.0);

    const double T = variable_array.temperature;

    const double S_v = 1 - S_L;
    const double D_vr = tortuosity_ * S_v;

    return base_diffusion_coefficient_ *
           std::pow(T / MaterialLib::PhysicalConstant::CelsiusZeroInKelvin,
                    exponent_) *
           D_vr;
}

PropertyDataType VapourDiffusionPMQ::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double S_L = std::clamp(variable_array.liquid_saturation, 0.0, 1.0);

    const double T = variable_array.temperature;

    if (variable == Variable::temperature)
    {
        const double S_v = 1 - S_L;
        const double D_vr = tortuosity_ * S_v;

        return exponent_ * base_diffusion_coefficient_ *
               std::pow(T / MaterialLib::PhysicalConstant::CelsiusZeroInKelvin,
                        exponent_ - 1.0) *
               D_vr / MaterialLib::PhysicalConstant::CelsiusZeroInKelvin;
    }

    if (variable == Variable::liquid_saturation)
    {
        return -base_diffusion_coefficient_ *
               std::pow(T / MaterialLib::PhysicalConstant::CelsiusZeroInKelvin,
                        exponent_) *
               tortuosity_;
    }

    OGS_FATAL(
        "VapourDiffusionPMQ::dValue is implemented for "
        "derivatives with respect to temperature or saturation only.");
}

}  // namespace MaterialPropertyLib
