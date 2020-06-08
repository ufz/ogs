/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
namespace PhreeqcIOData
{
struct Component;

std::vector<Component> createSolutionComponents(
    BaseLib::ConfigTree const& config, std::size_t const num_chemical_systems);
}  // namespace PhreeqcIOData
}  // namespace ChemistryLib
