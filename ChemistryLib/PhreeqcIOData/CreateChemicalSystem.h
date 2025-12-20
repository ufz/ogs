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
struct ChemicalSystem;

std::unique_ptr<ChemicalSystem> createChemicalSystem(
    BaseLib::ConfigTree const& config, MeshLib::Mesh& mesh);
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
