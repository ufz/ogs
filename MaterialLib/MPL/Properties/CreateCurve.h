// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <map>
#include <memory>
#include <string>

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{
class PiecewiseLinearInterpolation;
}

namespace MaterialPropertyLib
{
class Curve;

std::unique_ptr<Curve> createCurve(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);
}  // namespace MaterialPropertyLib
