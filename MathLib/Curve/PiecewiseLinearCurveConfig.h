// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

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
}  // namespace MathLib
