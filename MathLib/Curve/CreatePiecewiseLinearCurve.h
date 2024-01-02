/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Created on November 11, 2016, 10:49 AM
 */

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{
struct PiecewiseLinearCurveConfig
{
    std::vector<double> xs;
    std::vector<double> ys;
};

PiecewiseLinearCurveConfig parsePiecewiseLinearCurveConfig(
    BaseLib::ConfigTree const& config);

///  Create a curve
/// \param config   ConfigTree object has a tag of `<curve>`
template <typename CurveType>
std::unique_ptr<CurveType> createPiecewiseLinearCurve(
    BaseLib::ConfigTree const& config)
{
    auto [xs, ys] = parsePiecewiseLinearCurveConfig(config);
    return std::make_unique<CurveType>(std::move(xs), std::move(ys));
}
};  // namespace MathLib
