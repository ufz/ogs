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
class BishopsPowerLaw;
}

namespace MaterialPropertyLib
{
std::unique_ptr<BishopsPowerLaw> createBishopsPowerLaw(
    BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
