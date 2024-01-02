/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on March 4, 2021, 3:05 PM
 */

#include "LiquidViscosityVogels.h"

#include <cmath>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"

namespace MaterialPropertyLib
{

template <typename VogelsConstants>
PropertyDataType LiquidViscosityVogels<VogelsConstants>::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double T = variable_array.temperature;
    // Note: the constant of 1.e-3 is for the SI unit conversion.
    return 1.e-3 * std::exp(constants_.A + constants_.B / (constants_.C + T));
}

template <typename VogelsConstants>
PropertyDataType LiquidViscosityVogels<VogelsConstants>::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    if (variable != Variable::temperature)
    {
        OGS_FATAL(
            "LiquidViscosityVogels::dValue is implemented for "
            "derivatives with respect to temperature only.");
    }
    const double T = variable_array.temperature;
    const double f_buff = constants_.B / (constants_.C + T);
    // Note: the constant of 1.e-3 is for the SI unit conversion.
    return -1.e-3 * f_buff * std::exp(constants_.A + f_buff) /
           (constants_.C + T);
}

template class LiquidViscosityVogels<VogelsViscosityConstantsWater>;
template class LiquidViscosityVogels<VogelsViscosityConstantsCO2>;
template class LiquidViscosityVogels<VogelsViscosityConstantsCH4>;

}  // namespace MaterialPropertyLib
