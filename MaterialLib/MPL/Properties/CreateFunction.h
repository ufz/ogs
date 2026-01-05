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
class Function;
}

namespace MaterialPropertyLib
{
std::unique_ptr<Function> createFunction(
    BaseLib::ConfigTree const& config,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);
}  // namespace MaterialPropertyLib
