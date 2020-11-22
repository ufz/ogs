/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include <algorithm>
#include <cmath>

#include "KozenyCarmanModel.h"

namespace MaterialPropertyLib
{
PropertyDataType KozenyCarmanModel::value(
    MaterialPropertyLib::VariableArray const& variable_array,
    ParameterLib::SpatialPosition const& pos, double const t,
    double const /*dt*/) const
{
    double const phi = std::get<double>(variable_array[static_cast<int>(
        MaterialPropertyLib::Variable::porosity)]);
    auto const& k0 = _k0(t, pos);

    std::vector<double> k;
    k.reserve(k0.size());
    std::transform(k0.cbegin(), k0.cend(), std::back_inserter(k),
                   [&](auto const& k_component) {
                       return k_component *
                              std::pow((1 - _phi0(t, pos)[0]) / (1 - phi), 2) *
                              std::pow(phi / _phi0(t, pos)[0], 3);
                   });

    return fromVector(k);
}

}  // namespace MaterialPropertyLib
