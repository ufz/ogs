// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <optional>
#include <variant>
#include <vector>

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
namespace PhreeqcIOData
{
struct DensityBasedSurfaceSite;
struct MoleBasedSurfaceSite;

std::vector<std::variant<DensityBasedSurfaceSite, MoleBasedSurfaceSite>>
createSurface(std::optional<BaseLib::ConfigTree> const& config,
              MeshLib::Mesh& mesh);
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
