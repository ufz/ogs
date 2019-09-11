/**
 * \author Norbert Grunwald
 * \date   27.06.2018
 * \brief
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "MaterialLib/MPL/Properties/SaturationLiakopoulos.h"
#include "MaterialLib/MPL/Medium.h"
#include "MathLib/MathTools.h"

#include <algorithm>
#include <cmath>

namespace MaterialPropertyLib
{
/**
 */
PropertyDataType SaturationLiakopoulos::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& /*pos*/,
    double const /*t*/) const
{
    const double p_cap = std::get<double>(
        variable_array[static_cast<int>(Variable::capillary_pressure)]);

    if (p_cap < 0.)
        return 1.;

    return std::max(_residual_liquid_saturation,
                    1. - _parameter_a * std::pow(p_cap, _parameter_b));
}

PropertyDataType SaturationLiakopoulos::dValue(
    VariableArray const& variable_array, Variable const primary_variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/) const
{
    (void)primary_variable;
    assert((primary_variable == Variable::capillary_pressure) &&
           "SaturationLiakopoulos::dvalue is implemented for "
           " derivatives with respect to capillary pressure only.");

    const double p_cap = std::get<double>(
        variable_array[static_cast<int>(Variable::capillary_pressure)]);

    if ((p_cap < 0.) || (p_cap >= _p_cap_max))
    {
        return 0.;
    }
    return -_parameter_a * _parameter_b * std::pow(p_cap, _parameter_b - 1.);
}

PropertyDataType SaturationLiakopoulos::d2Value(
    VariableArray const& variable_array, Variable const primary_variable1,
    Variable const primary_variable2,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/) const
{
    (void)primary_variable1;
    (void)primary_variable2;
    assert((primary_variable1 == Variable::capillary_pressure) &&
           (primary_variable2 == Variable::capillary_pressure) &&
           "SaturationLiakopoulos::ddvalue is implemented for "
           " derivatives with respect to capillary pressure only.");

    const double p_cap = std::get<double>(
        variable_array[static_cast<int>(Variable::capillary_pressure)]);

    if ((p_cap < 0.) || (p_cap >= _p_cap_max))
    {
        return 0.;
    }

    return -_parameter_a * (_parameter_b - 1.) * _parameter_b *
           std::pow(p_cap, _parameter_b - 2.);
}

}  // namespace MaterialPropertyLib
