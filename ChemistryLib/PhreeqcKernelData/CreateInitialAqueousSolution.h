/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <string>
#include <vector>

namespace BaseLib
{
class ConfigTree;
}

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
class InitialAqueousSolution;

InitialAqueousSolution createInitialAqueousSolution(
    BaseLib::ConfigTree const& config,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map);
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
