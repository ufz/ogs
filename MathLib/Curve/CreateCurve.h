// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MathLib
{
class Curve;

/// Creates a curve from the `<curve>` subtree.
/// \param config ConfigTree object having the tag `<curve>`.
std::unique_ptr<Curve> createCurve(BaseLib::ConfigTree const& config);
}  // namespace MathLib
