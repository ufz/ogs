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
class PorosityFromMassBalance;
}

namespace ParameterLib
{
struct ParameterBase;
}

namespace MaterialPropertyLib
{
std::unique_ptr<PorosityFromMassBalance> createPorosityFromMassBalance(
    BaseLib::ConfigTree const& config,
    std::vector<std::unique_ptr<ParameterLib::ParameterBase>> const&
        parameters);
}  // namespace MaterialPropertyLib
