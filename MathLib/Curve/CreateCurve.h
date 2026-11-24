// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <string_view>

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{
class Curve;

/// Creates a curve from the `<curve>` subtree.
/// \param name Name of the curve, used in error messages.
/// \param config ConfigTree object having the tag `<curve>`.
std::unique_ptr<Curve> createCurve(std::string_view name,
                                   BaseLib::ConfigTree const& config);
}  // namespace MathLib
