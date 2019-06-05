/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "PhreeqcIO.h"

namespace BaseLib
{
class ConfigTree;
}

namespace ProcessLib
{
class Process;
}

namespace ChemistryLib
{
std::unique_ptr<PhreeqcIO> createPhreeqcIO(
    std::size_t const num_nodes,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map,
    BaseLib::ConfigTree const& config,
    std::string const& output_directory);
}  // namespace ChemistryLib
