/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
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
    const double S_L = std::clamp(
        std::get<double>(
            variable_array[static_cast<int>(Variable::liquid_saturation)]),
        0.0, 1.0);

    const double T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);
    const double phi =
        std::get<double>(variable_array[static_cast<int>(Variable::porosity)]);

    const double D_vr = tortuosity_ * phi * (1 - S_L);

    return 2.16e-5 *
           std::pow(T / MaterialLib::PhysicalConstant::CelsiusZeroInKelvin,
                    1.8) *
           D_vr;
}

PropertyDataType VapourDiffusionFEBEX::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double S_L = std::clamp(
        std::get<double>(
            variable_array[static_cast<int>(Variable::liquid_saturation)]),
        0.0, 1.0);

    const double T = std::get<double>(
        variable_array[static_cast<int>(Variable::temperature)]);
    const double phi =
        std::get<double>(variable_array[static_cast<int>(Variable::porosity)]);

    if (variable == Variable::temperature)
    {
        const double D_vr = tortuosity_ * phi * (1 - S_L);

        return 1.8 * 2.16e-5 *
               std::pow(T / MaterialLib::PhysicalConstant::CelsiusZeroInKelvin,
                        0.8) *
               D_vr / MaterialLib::PhysicalConstant::CelsiusZeroInKelvin;
    }
    if (variable == Variable::liquid_saturation)
    {
        return -2.16e-5 *
               std::pow(T / MaterialLib::PhysicalConstant::CelsiusZeroInKelvin,
                        1.8) *
               tortuosity_ * phi;
    }

    OGS_FATAL(
        "VapourDiffusionFEBEX::dValue is implemented for "
        "derivatives with respect to temperature or saturation only.");
}

}  // namespace MaterialPropertyLib
