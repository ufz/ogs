// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <optional>

#include "EquilibriumReactants.h"

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
std::unique_ptr<EquilibriumReactants> createEquilibriumReactants(
    std::optional<BaseLib::ConfigTree> const& config,
    MeshLib::Mesh const& mesh);
}
}  // namespace ChemistryLib
