/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "AqueousSolution.h"
#include "BaseLib/ConfigTree.h"
#include "CreateInitialAqueousSolution.h"
#include "InitialAqueousSolution.h"

namespace ChemistryLib
{
namespace PhreeqcKernelData
{
AqueousSolution createAqueousSolution(
    BaseLib::ConfigTree const& config,
    std::vector<std::pair<int, std::string>> const&
        process_id_to_component_name_map)
{
    //! \ogs_file_param{prj__chemical_system__solution__temperature}
    auto const temperature = config.getConfigParameter<double>("temperature");

    //! \ogs_file_param{prj__chemical_system__solution__pressure}
    auto const pressure = config.getConfigParameter<double>("pressure");

    //! \ogs_file_param{prj__chemical_system__solution__pe}
    auto const pe = config.getConfigParameter<double>("pe");

    auto const initial_aqueous_solution =
        createInitialAqueousSolution(config, process_id_to_component_name_map);

    return {temperature, pressure, pe, initial_aqueous_solution};
}
}  // namespace PhreeqcKernelData
}  // namespace ChemistryLib
