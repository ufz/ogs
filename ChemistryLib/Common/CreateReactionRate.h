/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <optional>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
template <typename ReactionRate>
std::vector<ReactionRate> createReactionRates(
    std::optional<BaseLib::ConfigTree> const& config);
}  // namespace ChemistryLib
