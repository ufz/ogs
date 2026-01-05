// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialPropertyLib
{
class Parameter;
}

namespace ParameterLib
{
struct ParameterBase;
}

namespace MaterialPropertyLib
{
std::unique_ptr<Parameter> createParameterProperty(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);
}  // namespace MaterialPropertyLib
