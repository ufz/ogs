/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "VermaPruessModel.h"

#include <algorithm>
#include <cmath>
#include <functional>

namespace MaterialPropertyLib
{
PropertyDataType VermaPruessModel::value(
    MaterialPropertyLib::VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    double const phi = std::get<double>(variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::porosity)]);

    // ratio of permeability to initial permeability
    auto const var_k = std::pow(std::max(0., phi - _phi_c(t, pos)[0]) /
                                    (_phi0(t, pos)[0] - _phi_c(t, pos)[0]),
                                _n(t, pos)[0]);

    auto const& k0 = _k0(t, pos);
    std::vector<double> k;
    k.reserve(k0.size());
    std::transform(
        k0.cbegin(), k0.cend(), std::back_inserter(k),
        std::bind(std::multiplies<double>(), std::placeholders::_1, var_k));

    return fromVector(k);
}

}  // namespace MaterialPropertyLib
