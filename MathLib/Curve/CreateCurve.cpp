// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include "CreateCurve.h"

#include "Curve.h"
#include "PiecewiseLinearCurveConfig.h"

namespace MathLib
{
std::unique_ptr<Curve> createCurve(BaseLib::ConfigTree const& config)
{
    auto [xs, ys] = parsePiecewiseLinearCurveConfig(config);
    return std::make_unique<Curve>(PiecewiseLinear{}, std::move(xs),
                                   std::move(ys));
}
}  // namespace MathLib
