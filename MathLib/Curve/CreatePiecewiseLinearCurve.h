// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string>
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

std::vector<double> readDoublesFromBinaryFile(const std::string& filename);

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
