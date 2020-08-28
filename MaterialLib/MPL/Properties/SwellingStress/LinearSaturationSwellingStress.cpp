/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on August 14, 2020, 8:19 AM
 */

#include "LinearSaturationSwellingStress.h"

#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/VariableType.h"

namespace MaterialPropertyLib
{
LinearSaturationSwellingStress::LinearSaturationSwellingStress(
    std::string name, double const coefficient,
    double const reference_saturation)
    : coefficient_(coefficient), reference_saturation_(reference_saturation)
{
    name_ = std::move(name);
}

PropertyDataType LinearSaturationSwellingStress::value(
    const VariableArray& variable_array,
    const ParameterLib::SpatialPosition& /*pos*/, const double /*t*/,
    const double dt) const
{
    // Sl <= S_max is guaranteed by the saturation property or
    // the saturation calculation.
    const double Sl = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    if (Sl < reference_saturation_)
    {
        return 0.0;
    }

    const double dS =
        dt *
        std::get<double>(
            variable_array[static_cast<int>(Variable::liquid_saturation_rate)]);

    return coefficient_ * dS;
}

PropertyDataType LinearSaturationSwellingStress::dValue(
    VariableArray const& variable_array, Variable const variable,
    ParameterLib::SpatialPosition const& /*pos*/, double const /*t*/,
    double const /*dt*/) const
{
    (void)variable;
    assert((variable == Variable::liquid_saturation) &&
           "LinearSaturationSwellingStress::dValue is implemented for "
           "derivatives with respect to liquid saturation only.");

    // Sl <= S_max is guaranteed by the saturation property or
    // the saturation calculation.
    const double Sl = std::get<double>(
        variable_array[static_cast<int>(Variable::liquid_saturation)]);

    return Sl < reference_saturation_ ? 0.0 : coefficient_;
}

}  // namespace MaterialPropertyLib
