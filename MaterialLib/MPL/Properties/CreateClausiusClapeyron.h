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
class ClausiusClapeyron;
}

namespace MaterialPropertyLib
{
std::unique_ptr<ClausiusClapeyron> createClausiusClapeyron(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
