// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ParameterLib
{
struct ParameterBase;
}

namespace MaterialPropertyLib
{
class Property;
}

namespace MaterialPropertyLib
{
std::unique_ptr<Property> createCubicLawPermeability(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);
}  // namespace MaterialPropertyLib
