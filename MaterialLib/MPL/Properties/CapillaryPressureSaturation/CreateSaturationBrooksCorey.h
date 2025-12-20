// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialPropertyLib
{
class SaturationBrooksCorey;
}

namespace MaterialPropertyLib
{
std::unique_ptr<SaturationBrooksCorey> createSaturationBrooksCorey(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
