// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <variant>

namespace MeshLib
{
class Mesh;
}

namespace BaseLib
{
class ConfigTree;
}

namespace NumLib
{
struct NoStabilization;
class IsotropicDiffusionStabilization;
class FullUpwind;
class FluxCorrectedTransport;

using NumericalStabilization =
    std::variant<NoStabilization, IsotropicDiffusionStabilization, FullUpwind,
                 FluxCorrectedTransport>;
}  // namespace NumLib

namespace NumLib
{
NumericalStabilization createNumericalStabilization(
    MeshLib::Mesh const& mesh, BaseLib::ConfigTree const& config);
}  // namespace NumLib
