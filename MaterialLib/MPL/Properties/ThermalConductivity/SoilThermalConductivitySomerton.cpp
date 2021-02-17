/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on February 17, 2021, 3:27 PM
 */

#include "SoilThermalConductivitySomerton.h"

#include <cmath>

#include "BaseLib/Error.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Utils/CheckVanGenuchtenExponentRange.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
PropertyDataType SoilThermalConductivitySomerton::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    const double S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    if (S_L <= 0.0)
    {
        return dry_thermal_conductivity_;
    }

    if (S_L > 1.0)
    {
        return wet_thermal_conductivity_;
    }

    return dry_thermal_conductivity_ +
           std::sqrt(S_L) *
               (wet_thermal_conductivity_ - dry_thermal_conductivity_);
}

PropertyDataType SoilThermalConductivitySomerton::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::liquid_saturation) &&
           "SoilThermalConductivitySomerton::dValue is implemented for "
           "derivatives with respect to liquid saturation only.");

    const double S_L = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    if (S_L <= 0.0 || S_L > 1.0)
    {
        return 0.0;
    }

    return 0.5 * (wet_thermal_conductivity_ - dry_thermal_conductivity_) /
           std::sqrt(S_L);
}

}  // namespace MaterialPropertyLib
