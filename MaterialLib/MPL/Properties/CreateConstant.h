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
class Constant;
}

namespace MaterialPropertyLib
{
std::unique_ptr<Constant> createConstant(BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
