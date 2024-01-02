/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

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
