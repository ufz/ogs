/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "TemperatureDependentDiffusion.h"

#include <algorithm>
#include <cmath>
#include <iterator>

#include "MaterialLib/PhysicalConstant.h"

namespace MaterialPropertyLib
{
void TemperatureDependentDiffusion::checkScale() const
{
    if (!std::holds_alternative<Component*>(scale_))
    {
        OGS_FATAL(
            "The property 'TemperatureDependentDiffusion' is "
            "implemented on the 'component' scale only.");
    }
}

PropertyDataType TemperatureDependentDiffusion::value(
    VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    auto const T = std::get<double>(variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::temperature)]);
    double const gas_constant = MaterialLib::PhysicalConstant::IdealGasConstant;
    double const Arrhenius_exponent =
        std::exp(Ea_ / gas_constant * (1 / T0_ - 1 / T));

    auto const D0_data = D0_(t, pos);
    std::vector<double> D;
    std::transform(D0_data.cbegin(), D0_data.cend(), std::back_inserter(D),
                   [&Arrhenius_exponent](double const D0_component)
                   { return D0_component * Arrhenius_exponent; });

    return fromVector(D);
}
}  // namespace MaterialPropertyLib
