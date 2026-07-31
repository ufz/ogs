// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

#include "Sigmoid.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MaterialPropertyLib
{
/// Create a Sigmoid property from XML configuration
std::unique_ptr<Sigmoid> createSigmoid(BaseLib::ConfigTree const& config);
}  // namespace MaterialPropertyLib
