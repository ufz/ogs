// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>

namespace BaseLib
{
class ConfigTree;
}

namespace MeshLib
{
class Mesh;
}  // namespace MeshLib

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct AqueousSolution;

std::unique_ptr<AqueousSolution> createAqueousSolution(
    BaseLib::ConfigTree const& config, MeshLib::Mesh& mesh);
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
