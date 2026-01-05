// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <memory>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
namespace SelfContainedSolverData
{
struct ChemicalReaction;

std::vector<std::unique_ptr<ChemicalReaction>> createChemicalReactionData(
    BaseLib::ConfigTree const& config);
}  // namespace SelfContainedSolverData
}  // namespace ChemistryLib
